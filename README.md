# Within-patient tumor–normal allelic expression in prostate cancer (CPGEA)

Analysis code for the manuscript *Paired tumor–normal analysis reveals selective strengthening of weak
allelic imbalance in prostate cancer*.

The cohort is CPGEA (Chinese Prostate Cancer Genome and Epigenome Atlas): matched tumor and
normal whole-genome sequencing and RNA-seq. The ASE atlas uses 118 tumors. The paired
tumor–normal analyses use 117 patients with a FACETS purity estimate and allelic coverage in
both tissues.

## Repository layout

| Directory | Manuscript section | Contents |
|---|---|---|
| `00_upstream/` | Methods | Allelic counting driver, FACETS, Beagle phasing, 16-bit count-overflow scan |
| `01_atlas_qc/` | Results 1, Fig. 1, Supp. Figs. S1–S3 | ASE calling QC, recurrence burden, expression dependence, GO enrichment, paired differential ASE (Supp. Table S16) |
| `02_attribution/` | Results 2, Fig. 2, Supp. Figs. S4–S6, S16 | SCNA, QTL, methylation and somatic-mutation attribution; event-level attribution |
| `03_paired_tn/` | Results 3 (H1), Fig. 3, Supp. Fig. S7 | Paired tumor–normal allelic deviation, mappability filter, robustness, MBASED |
| `04_tumor_state/` | Results 4 (H2–H4), Fig. 4A–E, Supp. Figs. S8, S10 | Δ-eQTL, baseline strata, epithelial composition and state, chromatin, cell-type interaction eQTL |
| `05_cis_by_trans/` | Results 5 (H5), Fig. 4F, Supp. Fig. S9, Supp. Table S17 | Transcriptome-wide cis-by-trans slope test, measurement-precision bias, TF sets, interaction-eQTL specificity |
| `06_kif18b/` | Results 6, Fig. 5, Supp. Fig. S11, Supp. Tables S8–S10, S15, S19, S20 | rs9914422 phasing, targeted screen, locus annotation, patient audit, external validation |
| `07_clinical/` | Supp. Figs. S13–S15, Supp. Tables S7, S11 | TCGA-PRAD, nested Cox, Gleason association |
| `08_figures/` | All figures | Publication figure scripts; each writes the statistics quoted in its legend |
| `config/` | | Environment, R packages, derived inputs, path remapping |

`MANIFEST.md` maps every script to the figure panel, table or statement it produces.

## Scope

This repository contains the code behind the submitted manuscript only. Exploratory analyses
that did not enter the manuscript, and earlier versions of analyses that were later redone,
are not included.

Several scripts carry a `_clean` suffix. They are the reported versions: they read the paired
site table after removal of 16-bit count-overflow records (`03_paired_tn/56_make_clean_sites.R`)
and cap depth weights at 100 reads. Where the uncleaned version is also included, it supplies a
sensitivity estimate quoted in a supplementary legend.

Code, comments and console messages are in English. Comments are kept to a minimum: scripts open
with a short header stating their purpose, inputs and outputs.

## Running the code

The scripts refer to inputs and outputs under the placeholder root `/path/to/data_root`, which
stands for the directory layout used when the reported results were produced (internal server
paths were replaced before release). To run the scripts elsewhere, copy the repository and point
the placeholder to your own data root:

```bash
cp -r pca-ase-code pca-ase-code-local
./config/remap_paths.sh pca-ase-code-local /your/data/root
```

Heavy steps (allelic counting, phasing, BayesPrism, permutation tests) were run as batch jobs
with 16–32 cores. Inputs that are produced once and reused are listed in `config/inputs.md`.

## Upstream method not redistributed

The per-gene phase-marginalised beta-binomial ASE model is `rna_allele_imbalance.R` from
Sen et al. (https://github.com/Arkosen/Allele-imbalance). The actual server caller is
`/path/to/data_root/user/github/Allele-imbalance-analysis-of-RNA-and-DNA/rna_allele_imbalance.R`
(SHA256 `c2f5237eeae518aff37d507d4cdfd553f8c1b75e19e98c8423f9dc8bfb5d1ccc`).
The archived `00_upstream/run_all_ASE_pipeline.sh` only aligns count-table labels and
integrates existing ASE results; it does not run STAR/WASP or the ASE test. Actual
STAR/WASP drivers were located under `/path/to/data_root/wasp/scripts/`
(`run_complete_pipeline_with_GQ_NEWDIR.sh` and the corresponding `_normal.sh`).
The historical
external directory is not a Git checkout, so its archived SHA256 identifies the
inspected version; a complete per-sample historical execution manifest is unavailable.

## Data access

CPGEA sequencing data are available through controlled access at the Genome Sequence Archive for
Human (https://ngdc.cncb.ac.cn/gsa-human/) under BioProject PRJCA001124 (GSA-Human study HRA000099). TCGA-PRAD data are available from
the Genomic Data Commons. No patient-level data are included in this repository.

## License

See `LICENSE`.
