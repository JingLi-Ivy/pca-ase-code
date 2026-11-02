#!/usr/bin/env python3
# 73 (2026-10-27): Supp. Table S21, cross-cohort candidate prioritization. Extracted verbatim from the v10 editorial
#     script that first wrote the table, so the selection rule is versioned with the analysis code.
# Inputs: saved CPGEA grade-expression and median-split recurrence results (gene_name, log2FC, pvalue, fdr, logrank_p,
#     cox_hr) and Supp. Table S7 (TCGA-PRAD, from 07_clinical/tcga_defB_final.R). Rule in each cohort: |log2FC| >= 0.5,
#     nominal P <= 0.05, log-rank P < 0.05 and concordant grade/recurrence direction; six genes pass in both.
# Usage: 73_s21_candidate_prioritization.py CPGEA_candidate_selection_input.csv TableS7_TCGA_prognosis.csv out.csv
import sys
import pandas as pd

cpgea, tcga, out = sys.argv[1:4]
cp = pd.read_csv(cpgea).set_index('gene_name')
tc = pd.read_csv(tcga).set_index('gene_name')
cp['selected_CPGEA'] = (cp.log2FC.abs() >= .5) & (cp.pvalue <= .05) & (cp.logrank_p < .05) & (cp.log2FC * (cp.cox_hr - 1) > 0)
tc['selected_TCGA'] = (tc.log2FC.abs() >= .5) & (tc.p_limma <= .05) & (tc.logrank_p < .05) & (tc.log2FC * (tc.HR_split - 1) > 0)
join = cp[['log2FC', 'pvalue', 'fdr', 'logrank_p', 'cox_hr', 'selected_CPGEA']].add_prefix('CPGEA_').join(
    tc[['log2FC', 'p_limma', 'fdr_limma', 'logrank_p', 'HR_split', 'fdr_logrank', 'selected_TCGA']].add_prefix('TCGA_'), how='outer')
join['selected_both'] = join.CPGEA_selected_CPGEA.fillna(False) & join.TCGA_selected_TCGA.fillna(False)
assert int(cp.selected_CPGEA.sum()) == 15 and int(tc.selected_TCGA.sum()) == 28
assert sorted(join[join.selected_both].index) == ['CNN1', 'DRD1', 'GLB1L3', 'KIF18B', 'NEK2', 'TOP2A']
join.reset_index().to_csv(out, index=False)
