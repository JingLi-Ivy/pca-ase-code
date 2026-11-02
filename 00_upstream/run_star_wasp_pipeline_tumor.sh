#!/usr/bin/env bash
# Per-sample allele-specific RNA read counting for a TUMOR sample (STAR + WASP).
#   A1-A2  STAR alignment (GRCh38/Ensembl index), sort and index
#   A3     per-chromosome heterozygous-SNP files from the sample's WGS genotypes
#   B      WASP mapping-bias filter: find_intersecting_snps.py -> STAR remap -> filter_remapped_reads.py
#   B1     MD tags, then WASP get_as_counts.py -> allele-specific counts
#   C      genotype error from GQ (10^(-GQ/10)) and SNP-to-gene (gene body) mapping
# Usage: run_star_wasp_pipeline_tumor.sh <SAMPLE_NAME>     (e.g. T32)
# Output: counts_for_R_pipeline_<SAMPLE>_with_GQ_geneBody.csv (input of the ASE caller)
# WASP (https://github.com/bmvdgeijn/WASP) and bedtools are external tools, not redistributed here.
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <SAMPLE_NAME>"
    echo "Example: $0 T32"
    exit 1
fi

source $HOME/miniconda3/bin/activate rnaseq_mutloh
export PATH="/path/to/data_root/user/github/bedtools2/bin:$PATH"

SAMPLE="$1"
VCF_SAMPLE="${SAMPLE}_WGS"
THREADS=12
BASE_DIR="/path/to/data_root/wasp"
cd "$BASE_DIR"

READ1="/path/to/data_root/raw/RNA_seq/tumor/${SAMPLE}_f1_trimmed.fq.gz"
READ2="/path/to/data_root/raw/RNA_seq/tumor/${SAMPLE}_r2_trimmed.fq.gz"
STAR_INDEX="/path/to/data_root/reference/hg38/STAR_index_ensembl_hg38"

STAR_DIR="${BASE_DIR}/star_ensembl_${SAMPLE}"
WASP_DIR="${BASE_DIR}/out_${SAMPLE}_ensembl"
SNP_DIR="${BASE_DIR}/snp_text_${SAMPLE}_het_nochr"
mkdir -p "$STAR_DIR" "$WASP_DIR" "$SNP_DIR"

export TMPDIR="/path/to/data_root/user/tmp"
mkdir -p "$TMPDIR"

log() { echo "[$(date '+%F %T')] $*"; }

if [ ! -f "$READ1" ]; then
    log "ERROR: input file not found: $READ1" >&2
    exit 1
fi
if [ ! -f "$READ2" ]; then
    log "ERROR: input file not found: $READ2" >&2
    exit 1
fi

log "Processing sample: $SAMPLE"
log "VCF sample name: $VCF_SAMPLE"

log "=== A1: STAR alignment ==="

UNSORTED_BAM="${STAR_DIR}/${SAMPLE}_Aligned.out.bam"
if [ ! -f "$UNSORTED_BAM" ]; then
    log "STAR alignment: ${SAMPLE} -> unsorted BAM"
    STAR \
        --runMode alignReads \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${STAR_DIR}/${SAMPLE}_" \
        --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD
    
    if [ -f "$UNSORTED_BAM" ]; then
        log "STAR alignment done: $UNSORTED_BAM"
    else
        log "ERROR: STAR output not found: $UNSORTED_BAM" >&2
        exit 1
    fi
else
    log "STAR output exists, skipping: $UNSORTED_BAM"
fi

log "=== A2: BAM sort and index ==="

SORTED_BAM="${STAR_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
if [ ! -f "$SORTED_BAM" ]; then
    log "samtools sort: $UNSORTED_BAM -> $SORTED_BAM"
    samtools sort -@ 12 -m 3G -T "$TMPDIR/${SAMPLE}.sort" -o "$SORTED_BAM" "$UNSORTED_BAM"
    log "samtools index: $SORTED_BAM"
    samtools index "$SORTED_BAM"
    log "BAM sorting done: $SORTED_BAM"
else
    log "Sorted BAM exists, skipping: $SORTED_BAM"
fi

log "=== A3: Building the SNP directory ==="

if [ ! -f "$SNP_DIR/1.snps.txt.gz" ]; then
    log "Creating SNP directory: $SNP_DIR"
    log "Extracting heterozygous SNPs of ${SAMPLE}..."
    bcftools view -s "$VCF_SAMPLE" -g het -v snps --threads 12 ${BASE_DIR}/combined.annotated.fixed.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | \
    awk 'BEGIN{OFS="\t"} { chr=$1; gsub(/^chr/,"",chr); if(length($3)==1 && length($4)==1){print chr,$2,$3,$4} }' | \
    awk 'BEGIN{OFS="\t"} { chr=$1; if(chr!=""){ print $2,$3,$4 > "'$SNP_DIR'/" chr ".snps.txt" } }'
    log "Compressing SNP files..."
    for file in "$SNP_DIR"/*.snps.txt; do
        [ -f "$file" ] && gzip "$file" && log "compressed: $file.gz"
    done
    log "Heterozygous SNPs of ${SAMPLE} extracted. Files: $(ls -1 "$SNP_DIR"/*.snps.txt.gz | wc -l)"
else
    log "SNP directory exists, skipping: $SNP_DIR"
fi

log "=== B: WASP remapping filter ==="

GENE_BODY_BED="${BASE_DIR}/gtf_gene_body_geneid.bed"
if [ -f "$GENE_BODY_BED" ]; then
    log "Reusing gene body BED: $GENE_BODY_BED"
else
    echo "[ERR] missing gene body BED: $GENE_BODY_BED" >&2; exit 1
fi
WASP_KEEP_BAM="${WASP_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.keep.bam"
WASP_FINAL_BAM="${WASP_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.keep.merged.sorted.bam"

if [ ! -f "$WASP_FINAL_BAM" ]; then
    PARALLEL_FIND=${PARALLEL_FIND:-1}
    FIND_N_PAR=${FIND_N_PAR:-4}

    if [ "$PARALLEL_FIND" = "1" ]; then
        log "Running find_intersecting_snps.py in parallel (jobs: $FIND_N_PAR)"
        CHR_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
        mkdir -p "$WASP_DIR/parallel"
        pjob=0
        part_keep_bams=()
        part_toremap_bams=()
        for chr in "${CHR_LIST[@]}"; do
            chr_bam="$TMPDIR/${SAMPLE}.chr${chr}.bam"
            out_chr_dir="$WASP_DIR/parallel/chr${chr}"
            mkdir -p "$out_chr_dir"
            samtools view -@ 12 -b "$SORTED_BAM" "$chr" -o "$chr_bam"
            samtools index "$chr_bam"
            (
                python /path/to/data_root/user/github/WASP/mapping/find_intersecting_snps.py \
                    --is_paired_end \
                    --snp_dir "$SNP_DIR" \
                    --output_dir "$out_chr_dir" \
                    --is_sorted "$chr_bam"
            ) &
            pjob=$((pjob+1))
            if [ $pjob -ge $FIND_N_PAR ]; then
                wait; pjob=0
            fi
        done
        wait
        while IFS= read -r -d '' kb; do part_keep_bams+=("$kb"); done < <(find "$WASP_DIR/parallel" -name "*.keep.bam" -print0)
        while IFS= read -r -d '' rb; do part_toremap_bams+=("$rb"); done < <(find "$WASP_DIR/parallel" -name "*.to.remap.bam" -print0)
        if [ ${#part_keep_bams[@]} -gt 0 ]; then
            log "Merging per-chromosome keep.bam -> $WASP_KEEP_BAM"
            samtools merge -@ 12 -f "$WASP_KEEP_BAM" "${part_keep_bams[@]}"
        else
            log "[WARN] no per-chromosome keep.bam found"
        fi
        MERGED_TOREMAP_BAM="$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.to.remap.bam"
        if [ ${#part_toremap_bams[@]} -gt 0 ]; then
            log "Merging per-chromosome to.remap.bam -> $MERGED_TOREMAP_BAM"
            samtools merge -@ 12 -f "$MERGED_TOREMAP_BAM" "${part_toremap_bams[@]}"
        else
            log "[INFO] per-chromosome to.remap is empty, creating an empty placeholder BAM"
            samtools view -H "$SORTED_BAM" | samtools view -b -o "$MERGED_TOREMAP_BAM"
        fi
    else
        log "WASP find_intersecting_snps.py (single job)"
        python /path/to/data_root/user/github/WASP/mapping/find_intersecting_snps.py \
            --is_paired_end \
            --snp_dir "$SNP_DIR" \
            --output_dir "$WASP_DIR" \
            --is_sorted "$SORTED_BAM"
        MERGED_TOREMAP_BAM="$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.to.remap.bam"
    fi

    if [ ! -f "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq1.gz" ]; then
        log "Merging per-chromosome fastq files..."
        cat "$WASP_DIR"/parallel/chr*/"${SAMPLE}".chr*.remap.fq1.gz > "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq1.gz"
        cat "$WASP_DIR"/parallel/chr*/"${SAMPLE}".chr*.remap.fq2.gz > "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq2.gz"
        
        if [ -s "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq1.gz" ]; then
            log "fastq merge succeeded, size: $(du -h "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq1.gz" | cut -f1)"
        else
            log "ERROR: fastq merge failed or file is empty" >&2
            exit 1
        fi
        log "fastq merge done"
    else
        log "Merged fastq files exist, skipping the merge"
    fi
    
    log "STAR remapping: ${SAMPLE}_Aligned.sortedByCoord.out.remap.fq1.gz + remap.fq2.gz"
    STAR \
        --runMode alignReads \
        --runThreadN 12 \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq1.gz" "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.fq2.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap." \
        --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD

    log "samtools sort: sorting remapped reads"
    samtools sort -@ 12 -m 2G -T "$TMPDIR/${SAMPLE}.remap.sort" -o "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.bam" "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.Aligned.out.bam"
    samtools index "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.bam"
    log "WASP filter_remapped_reads.py"
    python /path/to/data_root/user/github/WASP/mapping/filter_remapped_reads.py \
        "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.to.remap.bam" \
        "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.remap.bam" \
        "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.keep.remap.bam"

    log "WASP merge keep files"
    samtools merge -@ 12 -f "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.keep.merged.bam" \
        "$WASP_KEEP_BAM" "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.keep.remap.bam"
    samtools sort -@ 12 -m 2G -o "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.keep.merged.sorted.bam" \
        "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.keep.merged.bam"
    samtools index "$WASP_DIR/${SAMPLE}_Aligned.sortedByCoord.out.keep.merged.sorted.bam"
    log "WASP filtering done: $WASP_DIR"

    log "=== Read counts ==="
    echo "Original BAM reads: $(samtools view -c $SORTED_BAM)"
    echo "Keep BAM reads: $(samtools view -c $WASP_KEEP_BAM)"
    echo "Remap BAM reads: $(samtools view -c $MERGED_TOREMAP_BAM)"
    echo "Final BAM reads: $(samtools view -c $WASP_FINAL_BAM)"
else
    log "WASP output exists, skipping: $WASP_FINAL_BAM"
fi

log "=== B1: Adding MD tags and generating allele-specific counts ==="

REF_GENOME="/path/to/data_root/reference/ensembl_hg38.fa"
COUNTS_TXT="${BASE_DIR}/${SAMPLE}_allele_specific_counts.pipeline.txt"
COUNTS_NONZERO="${BASE_DIR}/${SAMPLE}_allele_specific_counts.pipeline.nonzero.txt"

if [ ! -f "$COUNTS_TXT" ]; then
    log "Adding MD tags to the BAM"
    samtools calmd -b "$WASP_FINAL_BAM" "$REF_GENOME" > "${WASP_FINAL_BAM%.bam}.with_md.unsorted.bam"
    samtools sort -@ 12 -m 3G -T "$TMPDIR/${SAMPLE}.withmd.sort" "${WASP_FINAL_BAM%.bam}.with_md.unsorted.bam" -o "${WASP_FINAL_BAM%.bam}.with_md.bam"
    samtools index "${WASP_FINAL_BAM%.bam}.with_md.bam"
    
    log "Running get_as_counts.py to generate allele-specific counts"
    python /path/to/data_root/user/github/WASP/mapping/get_as_counts.py --snp_dir "$SNP_DIR" "${WASP_FINAL_BAM%.bam}.with_md.bam" > "$COUNTS_TXT"
    
    log "Keeping non-zero counts"
    awk '($6+$7)>0' "$COUNTS_TXT" > "$COUNTS_NONZERO"
    
    log "Allele-specific counts done: $COUNTS_TXT"

    log "=== SNP coverage ==="
    echo "Total SNPs: $(zcat $SNP_DIR/*.snps.txt.gz | wc -l)"
    echo "SNPs with reads: $(awk '($6+$7)>0' $COUNTS_TXT | wc -l)"
    log "Non-zero count file: $COUNTS_NONZERO"
    log "Count file lines: $(wc -l < "$COUNTS_TXT" | tr -d ' ')"
    log "Non-zero count lines: $(wc -l < "$COUNTS_NONZERO" | tr -d ' ')"
else
    log "Allele-specific count file exists, skipping: $COUNTS_TXT"
    if [ ! -f "$COUNTS_NONZERO" ]; then
        log "Generating the non-zero count file"
        awk '($6+$7)>0' "$COUNTS_TXT" > "$COUNTS_NONZERO"
    fi
fi

log "=== C: GQ calculation and final output ==="

VCF_HET="${BASE_DIR}/${SAMPLE}_het_pipeline.vcf.gz"
if [ ! -f "$VCF_HET" ]; then
    log "bcftools: extracting heterozygous SNPs: $VCF_HET"
    bcftools view -s "$VCF_SAMPLE" -g het -v snps --threads 12 -Oz -o "$VCF_HET" ${BASE_DIR}/combined.annotated.fixed.vcf.gz
    tabix -p vcf "$VCF_HET"
else
    log "Heterozygous SNP VCF exists: $VCF_HET"
fi

GQ_MAP="${BASE_DIR}/snp_gq_mapping_${SAMPLE}.txt"
log "Writing GQ map: $GQ_MAP"
bcftools query -s "$VCF_SAMPLE" -f '%CHROM\t%POS\t[%GQ]\n' "$VCF_HET" \
| awk 'BEGIN{OFS="\t"}{c=$1; sub(/^chr/,"",c); g=$3; if(g==""||g==".") g=1; print c ":" $2, g}' \
> "$GQ_MAP"
log "GQ_MAP lines: $(wc -l < "$GQ_MAP" | tr -d ' ')"

log "=== Genotype-quality distribution of heterozygous SNPs ==="
bcftools query -s "$VCF_SAMPLE" -f '[%GQ]\n' "$VCF_HET" | \
    awk '{if($1>=20) high++; else low++} END{print "GQ>=20:", high, "GQ<20:", low}'

COUNTS_NONZERO="${BASE_DIR}/${SAMPLE}_allele_specific_counts.pipeline.nonzero.txt"
CORE_OUT="${BASE_DIR}/counts_core_pipeline_${SAMPLE}_with_GQ.csv"
log "Recomputing core table: $CORE_OUT"
awk -F'[ \t]+' -v OFS=',' -v map="$GQ_MAP" '
BEGIN{ while((getline l < map)>0){ split(l,a,"\t"); g[a[1]]=a[2] } close(map) }
{
  key=$1 ":" $2; pos=$2; ref=$6; alt=$7; oth=$8;
  ge=(key in g)?g[key]:1; if(ge<=0) ge=1;
  geno=10^(-ge/10);
  printf "%s,%s,%s,%s,%.6e,%s\n", key,pos,ref,ref+alt,geno,oth;
}' "$COUNTS_NONZERO" > "$CORE_OUT"

GENE_BODY_BED="${BASE_DIR}/gtf_gene_body_geneid.bed"
if [ ! -f "${GENE_BODY_BED}" ]; then
  echo "[ERR] missing gene body BED: ${GENE_BODY_BED}" >&2
  echo "Generate it first with 00_upstream/build_gene_body_bed.sh" >&2
  exit 1
fi

SNPS_BED="${BASE_DIR}/${SAMPLE}_snps_pipeline.bed"
awk -F'[ \t]+' 'BEGIN{OFS="\t"} {c=$1; sub(/^chr/,"",c); pos=$2; print c, pos-1, pos, c ":" pos}' \
  "${COUNTS_NONZERO}" | sort -k1,1V -k2,2n -u > "${SNPS_BED}"

MAPPING_BODY="${BASE_DIR}/snp_to_gene_pipeline_${SAMPLE}_geneBody.csv"
bedtools intersect -a "${SNPS_BED}" -b "${GENE_BODY_BED}" -wa -wb \
| awk -v OFS=',' '{print $4, $8}' \
| sort -u > "${MAPPING_BODY}"

FINAL_OUT_BODY="${BASE_DIR}/counts_for_R_pipeline_${SAMPLE}_with_GQ_geneBody.csv"
{
  echo 'gene_id,start,ref.matches,N,genotype.error,errors';
  awk -F',' -v OFS=',' 'NR==FNR{
      m[$1] = (m[$1] ? m[$1] ";" $2 : $2); 
      next
    }{
      key = $1; 
      if(!(key in m)) next; 
      n = split(m[key], g, ";");
      for(i=1;i<=n;i++) print g[i], $2, $3, $4, $5, $6;
  }' "${MAPPING_BODY}" "${CORE_OUT}";
} > "${FINAL_OUT_BODY}"

log "=== Done ==="
log "Sample: $SAMPLE"
log "Final output: ${FINAL_OUT_BODY}"
log "Lines: $(wc -l < "${FINAL_OUT_BODY}" | tr -d ' ')"
echo "Done: ${FINAL_OUT_BODY}"
