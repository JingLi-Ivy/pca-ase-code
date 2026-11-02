#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l cpu=32:mem=256G
#$ -q fast
#$ -N facets
#$ -o /path/to/data_root/variants/batch_43/facets/logs/facets.out
#$ -e /path/to/data_root/variants/batch_43/facets/logs/facets.err
#$ -V

echo "============ Debug info ============"
echo "Host: $(hostname)"
echo "User: $(whoami)"
echo "Directory: $(pwd)"
echo "=================================="

if [ $# -lt 3 ]; then
    echo "Error: missing required arguments"
    echo "Usage: $0 <patient_id> <tumor_bam> <normal_bam>"
    exit 1
fi

PATIENT_ID="$1"
TUMOR_BAM="$2"
NORMAL_BAM="$3"

echo "Patient ID: $PATIENT_ID"
echo "Tumor BAM: $TUMOR_BAM"
echo "Normal BAM: $NORMAL_BAM"

FACETS_DIR="/path/to/data_root/variants/batch_43/facets"
PILEUP_DIR="$FACETS_DIR/pileup"
RESULTS_DIR="$FACETS_DIR/results"
PILEUP_FILE="$PILEUP_DIR/${PATIENT_ID}.pileup.csv.gz"
OUTPUT_PREFIX="$RESULTS_DIR/${PATIENT_ID}"

mkdir -p "$PILEUP_DIR"
mkdir -p "$RESULTS_DIR"

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_mutloh

if [ ! -f "$TUMOR_BAM" ]; then
    echo "Error: tumor BAM file not found: $TUMOR_BAM"
    exit 1
fi

if [ ! -f "$NORMAL_BAM" ]; then
    echo "Error: normal BAM file not found: $NORMAL_BAM"
    exit 1
fi

if [ ! -f "$PILEUP_FILE" ]; then
    echo "Generating pileup file: $PILEUP_FILE"
    snp-pileup -g -q15 -Q20 -P100 -r25,0 \
        $HOME/miniconda3/envs/rnaseq_mutloh/lib/R/library/facets/extdata/dbsnp_137.hg19.vcf.gz \
        "$PILEUP_FILE" "$NORMAL_BAM" "$TUMOR_BAM"
else
    echo "pileup file already exists, reusing: $PILEUP_FILE"
fi

echo "Starting FACETS analysis..."
Rscript "$FACETS_DIR/run_facets.R" "$PILEUP_FILE" "$OUTPUT_PREFIX"

echo "FACETS analysis done: $PATIENT_ID"
