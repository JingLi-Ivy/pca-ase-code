# Environment

## R

R 4.4.1. Packages are listed in `r_packages.txt`. Bioconductor packages were installed with
`BiocManager::install()`; `facets` was installed from its GitHub repository.

Analyses were run inside a conda environment named `r_env`:

```bash
conda activate r_env
Rscript <script>.R
```

## Python

Python 3 is called only by the upstream allelic-counting driver (`00_upstream/run_all_ASE_pipeline.sh`).

## External tools

| Tool | Used for |
|---|---|
| STAR | RNA-seq alignment |
| WASP | Removal of reference-mapping bias before allelic counting |
| GATK (HaplotypeCaller, CalculateGenotypePosteriors) | Germline genotyping |
| Mutect2 | Paired tumour/normal somatic variant calling |
| bcftools, vcftools, SnpSift, PLINK | Variant processing, filtering, LD |
| FACETS | Allele-specific copy number, purity, ploidy |
| ANNOVAR | Annotation of panel-of-normals-filtered somatic variants |
| Beagle 5.5 (1000 Genomes GRCh38 reference panel) | Statistical phasing of germline SNPs |
| Bismark | WGBS alignment and methylation calling |
| samtools | BAM handling |

## Compute

Resource-intensive steps were submitted to a Sun Grid Engine cluster, for example:

```bash
qsub -cwd -l cpu=32:mem=256G --running_timeout=720h -V -q fast -N <job_name> <script>.sh
```

Scripts that are wrapped in a `.sh` submission file assume this scheduler; they run unchanged
under `Rscript` on a machine with sufficient memory. The heaviest steps (the genome-wide
cis–trans screen and the event-level attribution) require roughly 128–256 GB of RAM.
