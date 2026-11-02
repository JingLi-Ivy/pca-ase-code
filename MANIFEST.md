# Manifest: script → manuscript item

Figure and table numbers refer to the submitted manuscript. "Supp." = Supplemental.
Scripts are listed in the order they run within each directory.

## 00_upstream — allelic counts, copy number, phasing

| Script | Produces |
|---|---|
| `run_star_wasp_pipeline_tumor.sh`, `run_star_wasp_pipeline_normal.sh` | Per-sample STAR alignment, WASP mapping-bias filtering and allele-specific read counting (tumor and adjacent normal), with genotype error from GQ and SNP-to-gene mapping. Produces the allelic count tables used by every ASE analysis (Methods). |
| `run_all_ASE_pipeline.sh` | Count-table label alignment and integration of existing ASE results; it does not run STAR/WASP or the ASE test (see README). |
| `run_facets.sh`, `run_facets.R` | FACETS allele-specific copy number, purity and ploidy (hg38). Defines balanced-diploid sites and SCNA states used throughout. |
| `08_kif18b_phasing.sh`, `08a_align_alleles.R` | Beagle 5.5 phasing of the KIF18B locus against 1000 Genomes GRCh38, after aligning CPGEA REF/ALT to the reference panel. Input to `06_kif18b/08b_phased_allelic_ratio.R`. |
| `31_phase_genome.sh`, `31a_align_fast.R` | Genome-wide Beagle 5.5 phasing of germline SNPs. Input to `03_paired_tn/37d_mbased_phased.R`. |
| `55_count_overflow_scan.sh`, `55b_overflow_genes.R` | Scan of all allele-count files for 16-bit integer overflow (depth > 32,767) and the genes at risk; none is a PCa-ASE gene (Methods). |

## 01_atlas_qc — ASE atlas (Results 1, Fig. 1)

| Script | Produces |
|---|---|
| `00_pca_ase_catalogue.R` | Supp. Table S2 and the 2,378-gene PCa-ASE catalogue used throughout (definition in Methods; Sen et al. 2022 criteria applied per gene symbol). |
| `07_ASE_QC_panel.R` | ASE calling QC: depth distribution, reference-allele fraction, beta-binomial null calibration, depth-threshold sensitivity. |
| `61a_supp_qc_data.R` | Plot data for Supp. Fig. S1 A–C, recomputed with the procedure and seed of `07_ASE_QC_panel.R`. |
| `29_fig1d_threeset_burden.R` | Fig. 1D. Mean ASE magnitude per gene, tumor versus normal, for all genes, Cancer Gene Census genes and imprinted controls. |
| `53_s4c_expr_strat.R` | Supp. Fig. S2B. Frequency of significant ASE by expression tertile. |
| `68_go_recheck_allterms_fdr.R` | Supp. Fig. S3 data. topGO enrichment of the 2,378 PCa-ASE genes against the background of genes with heterozygous-SNP coverage; Benjamini–Hochberg over all tested terms. |
| `paired_diff_ase.R` | Supp. Table S16, columns without coverage control. |
| `paired_diff_ase_depthmatched_gw.R`, `paired_diff_ase_dm_analyze.R` | Supp. Table S16, coverage-matched columns: shared heterozygous sites, hypergeometric downsampling to equal depth in both tissues, per-gene paired tests. |

## 02_attribution — source attribution (Results 2, Fig. 2)

| Script | Produces |
|---|---|
| `04_Fig2_ASE_CNV_analysis.R` | Fig. 2B, E, F and Supp. Table S3. ASE–SCNA Spearman correlation for 1,796 genes; ATP2C1 and MCM9 examples. |
| `fig2C.R` | Fig. 2C. Deletion- versus amplification-associated genes. |
| `rerun_S4_KEGG_in_RStudio.R` | Supp. Table S4. KEGG enrichment of deletion-driven ASE genes against the tested-gene background. |
| `integrate_pca_ase_mechanisms.R` | Mutually exclusive evidence combinations for the 2,378 PCa-ASE genes (former pie chart; not in the submitted version, where the UpSet plot in Fig. 2D shows the same counts). |
| `Pca_ase_mechanism_integration_v2.R` | Fig. 2D and Supp. Table S5. Non-exclusive source annotation (UpSet). |
| `06_qtl_ld_attribution.R` | QTL attribution sensitivity quoted in Results 2: LD (r² ≥ 0.8) instead of positional identity, and gene-level QTL status. |
| `eqtl_concordance_full.R` | Supp. Fig. S4 C. eQTL–ASE direction concordance. |
| `66_s6df_ase_qtl_enrichment.R` | Supp. Fig. S4 A–B. MAF-matched and expression-stratified QTL enrichment among ASE SNPs. |
| `08_methylation_triangulation.R` | Supp. Fig. S5 and Supp. Table S13. Promoter plus gene-body DMR methylation–ASE–expression triangulation; current v10 event intersection. |
| `62d_pon_exonic_extract.sh`, `62b_s6_mutation_data.R`, `62e_s6_pon_summary.R` | Supp. Fig. S6 A, C. Panel-of-normals-filtered somatic variants in the 118 ASE-cohort tumors; consequence counts and overlap with SCNA-associated genes. |
| `62f_s6b_pon_permutation.R` | Supp. Fig. S6 B. Missense enrichment among copy-neutral ASE cells, 1,000 within-gene permutations. |
| `72_DNA_RNA_harmonized.R` | Corrected v10 DNA–RNA concordance: explicit REF/ALT identity, unique allele-site records and patient bootstrap. Writes isolated revalidation outputs. |
| `02_event_attribution_domseg.R` | Supp. Fig. S16 and Supp. Table S12. Event-level attribution of 14,622 significant tumor ASE events by FACETS state, with confidence tiers. |

## 03_paired_tn — paired tumor–normal comparison, H1 (Results 3, Fig. 3)

| Script | Produces |
|---|---|
| `18_delta_ase_pilot.R` | Paired site table: germline heterozygous SNPs, tumor balanced diploid, ≥ 10 reads in both tissues. |
| `19_delta_ase_full.R` | Leave-one-out orientation and the per-SNP paired test framework used by all downstream scripts. |
| `56_make_clean_sites.R` | Removes 16-bit count-overflow records and overflow-risk genes from the paired site table. All `_clean` scripts read its output. |
| `27_delta_ase_mapfilter_clean.R` | **H1 as reported.** Mappability filter (segmental duplications, ENCODE blacklist, MHC, pseudogenes, consistent near-monoallelic sites); patient-level tumor/normal ratio with patient bootstrap (Fig. 3A); per-SNP quasi-binomial test with 20 within-patient label swaps (Fig. 3B, Supp. Table S14); split-sample direction (Fig. 3C). |
| `27_delta_ase_mapfilter.R` | The same analysis before overflow cleaning; the original-weights estimate in Supp. Fig. S7A. |
| `27u_h1_unfiltered_clean.R` | H1 with and without the mappability filter. |
| `27x_overflow_check.R`, `27y_weight_sensitivity.R` | Supp. Fig. S7A. Sensitivity to overflow removal and to depth-weight caps (none, 1,000, 100, 30 reads, unweighted). |
| `20_delta_ase_dna_check.R` | Fig. 3D. Tumor versus normal DNA allelic deviation at tumor-stronger SNPs (copy-number check). |
| `54_tn_identity_check.R` | Supp. Fig. S7B. Fraction of near-monoallelic sites per patient in tumor and normal (sample identity). |
| `37b_mbased_rerun.R` | MBASED with pseudo-phasing (sensitivity, Supp. Fig. S7). |
| `37d_mbased_phased.R` | Fig. 3E, F and Supp. Fig. S7C. MBASED with Beagle phase; tumor versus normal major-allele frequency per patient. |

## 04_tumor_state — H2–H4 (Results 4, Fig. 4)

| Script | Produces |
|---|---|
| `32_bayesprism_v2_counts.R` | BayesPrism deconvolution of tumor and normal RNA-seq with the GSE185344 single-cell reference (cell-type and epithelial-state fractions). Run with `DECON_OUT=.../decon_v2d` for the state-resolved fractions. |
| `34c_delta_eqtl_dedup_lmm_clean.R` | Fig. 4A, H2. Within-patient Δ-eQTL at GTEx v8 prostate lead eQTLs, read de-duplication, linear mixed model. |
| `38_delta_decompose_clean.R` | Fig. 4C. Tumor/normal ratio by GTEx v8 prostate eGene status (and expression tertile). |
| `46_h3_trend_clean.R` | Fig. 4B, H3. Tumor/normal ratio by normal-tissue baseline tertile over 20 random split-sample partitions. |
| `51_h4_splits_clean.R` | Fig. 4D, H4. Cross-half orientation by the tumor- or normal-preferred allele, same split design. |
| `40_epithelial_dilution_clean.R` | Fig. 4E (all patients, epithelium-rich normal, similar epithelial fraction). Epithelial dilution model. |
| `41_epithelial_state_clean.R` | Fig. 4E (similar basal share, similar tumor-state share). Within-epithelium state composition. |
| `42_awakened_chromatin_tf.R`, `43_genelevel_chromatin.R` | Results 4. Weak-baseline gain versus TCGA PRAD ATAC-seq peaks (gene-level, length-adjusted) and allele-specific binding. |
| `48_celltype_ieqtl.R` | Supp. Fig. S10. Genotype × cell-type-fraction interaction eQTL at GTEx v8 prostate leads, permutation null. |

## 05_cis_by_trans — H5 (Results 5, Fig. 4F)

| Script | Produces |
|---|---|
| `feasibility.R` | Per-gene, per-sample aggregate of the Sen statistic and allelic depth in tumors. |
| `prep_normal.R` | The same aggregate for normal tissue. |
| `26e_tn_cistrans_slope_ndepth.R` | **H5 as reported.** Tumor-minus-normal cis-by-trans slope difference (heterozygous versus homozygous), covariate models, 1,000 whole-patient permutations (Fig. 4F; Supp. Fig. S9; Supp. Table S17 A). |
| `49_reliability_bias.R` | Supp. Fig. S9 A. Measurement-precision bias of the correlation-difference statistic, simulated on the real data structure. |
| `50_h5_bias_expectation.R` | Expected H5 under no cis-by-trans given that bias (Supp. Table S17). |
| `52_prostate_tf_set.R` | Prostate lineage and driver TF sets versus expression-matched TF sets (Supp. Table S17). |
| `28d_ieqtl_specificity_1000draws.R` | Historical CPGEA-lead specificity analysis; not the source of final v10 Supp. Table S17F. |
| `28f_ieqtl_specificity_gtex_paired.R` | **As reported in Results 5, Supp. Fig. S9D and Supp. Table S17F.** The test at independent GTEx v8 prostate lead eQTLs in the paired patients. |

## 06_kif18b — the KIF18B locus (Results 6, Fig. 5)

| Script | Produces |
|---|---|
| `08b_phased_allelic_ratio.R` | Phased, tumor-DNA-corrected C-haplotype fraction of KIF18B RNA in rs9914422 heterozygotes (Supp. Table S10). Plotted by `08_figures/25_kif18b_phase_panel.R`. |
| `11_targeted_screen_v8r2.R` | Supp. Tables S9 and S15. Targeted variant × regulator interaction screen on target-gene expression (133 tumors); all regulators at rs9914422–KIF18B. |
| `FOXM1_ASE_analysis.R`, `40_fig5c_clean.R` | Supp. Table S8. Significant ASE frequency of FOXM1 target genes by FOXM1 expression tertile. |
| `22_kif18b_locus_figure.R`, `adastra_extract.R` | Locus annotation inputs for Fig. 5A (gene models, ENCODE cCREs, reporter fragment) and ADASTRA v6.1 allele-specific binding at rs9914422. |
| `FigS5_rs9914422_LD.R` | Supp. Fig. S11A. Pairwise LD at the KIF18B locus from CPGEA germline genotypes. |
| `04_kif18b_patient_audit.R` | Supp. Table S20. Patient-level KIF18B evidence. |
| `53_kif18b_ase_gleason_matched.R` | KIF18B ASE versus Gleason against expression-matched genes (Supp. Fig. S11). |
| `21_kif18b_beyond_proliferation.R` | Supp. Fig. S13B. KIF18B hazard ratio after adjustment for proliferation scores. |
| `10b_external_validation.R` | Supp. Table S19. Pre-specified external validation in four cohorts. |

## 07_clinical — clinical analyses (Supp. Figs. S13–S15)

| Script | Produces |
|---|---|
| `candidates_tcga_batch.R` | TCGA-PRAD expression matrix and clinical table (workspace for the next script). |
| `tcga_defB_final.R` | Supp. Table S7 and Supp. Fig. S13A. TCGA-PRAD limma by grade, log-rank and Cox for biochemical recurrence. |
| `09_nested_risk_model.R`, `31_nested_risk_model_pub.R` | Supp. Table S11 and Supp. Fig. S15A–B. Nested Cox models with optimism-corrected C-index. |
| `ase_clinical_robust.R`, `ase_gleason_true_copyneutral.R` | Supp. Fig. S14. Global allelic-imbalance magnitude versus Gleason grade, all sites and balanced-diploid regions. |
| `35_delta_clinical_clean.R` | Patient-level Δ-ASE and Δ-eQTL versus clinical features (exploratory; Supplemental Results). Its patient table supplies purity for Supp. Fig. S15C–D. |
| `45c_awakened_gain_clinical.R` | Supp. Fig. S15 C–D. Pre-specified clinical test of the weak-baseline gain. |
| `73_s21_candidate_prioritization.py` | Supp. Table S21. Cross-cohort candidate prioritization (CPGEA saved results joined with Supp. Table S7); reproduces the submitted table byte for byte. |
| `11_celltype_localization.R`, `12_reannotate_localization.R` | KIF18B, FOXM1, NEK2 and TOP2A across cell types in six single-cell datasets. Not included in the submitted version. |

## Supplemental tables without a row above

| Table | Source |
|---|---|
| S1 | Clinical characteristics compiled from the CPGEA clinical annotation (`CPGEA_Clinical_Info.tsv`); descriptive, no analysis script. |
| S2 | `01_atlas_qc/00_pca_ase_catalogue.R`. Reconstructs the 2,378-gene PCa-ASE catalogue from `ase_data.RData` (saved 2025-12-19) and reproduces the submitted table exactly (all columns and ranks). |
| S13B | `02_attribution/08_methylation_triangulation.R`, intersection with the event table of `02_attribution/02_event_attribution_domseg.R`. |
| S16B | `01_atlas_qc/paired_diff_ase_depthmatched_gw.R` (`gw_depthmatched_pergene_D8.csv`, complete depth-8 testing family). |
| S18 | siRNA and qRT-PCR primer sequences from the wet-lab records; no analysis script. |
| S22 | `02_attribution/72_DNA_RNA_harmonized.R` (`DNA_RNA_allele_harmonized_summary.csv`). |

## 08_figures — publication figures

Each script reads saved results only and writes the numbers quoted in its legend to a
`*_legend_stats.txt` file next to the figure.

Supplemental figures were renumbered after these scripts were written. Output file names keep the
working numbers: `SuppFig_S15`, `SuppFig_S16` and `SuppFig_S17` are submitted Supplemental Figs. S14,
S15 and S16; `SuppFig_S14` (single-cell localization) and the `suppS18` GO file name predate the final
numbering (the GO dot plot is Supplemental Fig. S3). The tables above use the submitted numbers.

| Script | Produces |
|---|---|
| `fig_theme.R` | Shared style (Arial 8 pt, final print size). |
| `58_fig1_pub.R` | Fig. 1 C–F. |
| `59a_fig2A_panel.R` | Fig. 2A. FACETS copy-number landscape of the 118 ASE-cohort tumors (5-Mb bins, ordered by fraction of genome altered, log ratio capped at ±1). Panels B–G of Fig. 2 come from the `02_attribution` scripts. |
| `57_fig3_fig4_pub.R` | Fig. 3 A–F and Fig. 4 A–F. |
| `25_kif18b_phase_panel.R` | Fig. 5B. Pooled beta-binomial C fraction with DNA offset. |
| `60_fig5_pub.R` | Fig. 5C (dual-luciferase reporter; two-way ANOVA on log activity) and Fig. 5E (siFOXM1 qRT-PCR; the one-way ANOVA P values in the legend were computed in Excel). |
| `70_fig5d_suppS17.R` | Fig. 5D (FOXM1 → KIF18B in 118 tumors) and Supp. Fig. S16. |
| `72_fig5_assemble.R` | Fig. 5 assembled; draws panel A (locus, browser-track style) and panel F (model summary). |
| `61_supp_figs_A.R` | Supp. Figs. S1, S2, S4, S5, S14. |
| `69_suppS18_go_dotplot.R` | Supp. Fig. S3 (GO). |
| `64_supp_S6.R` | Supp. Fig. S6. |
| `62_supp_figs_B.R` | Supp. Figs. S7–S10. |
| `63_supp_figs_C.R` | Supp. Figs. S11–S13, S15 (and the single-cell panel not in the submitted version). |
| `67_supp_tables_S14_S17.R` | Supp. Tables S14 and S17, with built-in consistency assertions; S17F uses v10 28f. |
| `71_export_v10_S6.R` | Corrected Supp. Table S6 from the v10 PoN/ANNOVAR input, with sample-gene enrichment categories; full input companion export. |

Fig. 1A–B are schematics drawn by hand.

## 2026-10-13 figure layout pass (v10 submission figures)

These scripts supersede the corresponding figure scripts above for the submitted figures. They change layout, labels and colours only; the plotted data and every statistic quoted in the legends are unchanged (checked against the previous `*_legend_stats.txt`).

| Script | Produces |
|---|---|
| `08_figures/58b_fig1_pub.R` | Fig. 1 C–F (A, B are hand-drawn schematics) |
| `08_figures/59b_fig2_pub.R` | Fig. 2 A–F (A from `59a`; UpSet replaces the former pie chart, so former E–G are now D–F) |
| `08_figures/57b_fig3_fig4_pub.R` | Fig. 3 and Fig. 4 |
| `08_figures/61b_supp_figs_A.R` | Supp. Figs. S1, S2, S4, S5, S14 (S3 is `69_suppS18_go_dotplot.R`) |
| `08_figures/62b_supp_figs_B.R` | Supp. Figs. S7–S10 |
| `08_figures/63b_supp_figs_C.R` | Supp. Figs. S11–S13, S15 (and the single-cell panel not in the submitted version) |
