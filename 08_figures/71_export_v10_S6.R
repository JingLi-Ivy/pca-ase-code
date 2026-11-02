#!/usr/bin/env Rscript
# v10 correction: export the PoN/ANNOVAR variants underlying 62d/62e/62f, not the legacy VEP table.
suppressPackageStartupMessages(library(data.table))
V <- '/path/to/data_root/ASE/code/v8_round2'
P <- '/path/to/data_root/ASE'
O <- file.path(V, 'supp_data')
x <- fread(file.path(O, 'S6_pon_exonic_variants.tsv'))
x[, gene_ANNOVAR_original := gene]
x[, gene := sub('(\\\\x3b|[,;]).*', '', gene)]
g <- fread(file.path(P, 'code/PCa_ASE_genes.csv'))[, .(gene_symbol, g = sub('\\..*', '', gene_id))]
x <- merge(x, unique(g), by.x = 'gene', by.y = 'gene_symbol', all.x = TRUE)
x[, in_PCa_ASE_catalogue := !is.na(g)]
cells <- fread(file.path(O, 'S6B_pon_per_cell.csv'))[, .(sample, g, category)]
stopifnot(!anyDuplicated(cells, by = c('sample', 'g')))
x <- merge(x, cells, by = c('sample', 'g'), all.x = TRUE)
x[is.na(category), category := fifelse(in_PCa_ASE_catalogue, 'not_evaluable_for_enrichment', 'outside_PCa_ASE_catalogue')]
x[, variant_source := 'PoN-filtered ANNOVAR; v10 62d/62e/62f']
stopifnot(nrow(x) == 12984L, uniqueN(x$sample) == 118L)
y <- x[in_PCa_ASE_catalogue == TRUE]
stopifnot(nrow(unique(y[exonic_func == 'nonsynonymous_SNV', .(sample, gene)])) == 520L)
fwrite(y, file.path(O, 'TableS6_somatic_mutations.csv'))
fwrite(x, file.path(O, 'TableS6_all_genes_PoN_input.csv'))
cat('S6 variant rows:', nrow(y), '; tumors with PCa-ASE gene variants:', uniqueN(y$sample), '; assayed cohort:', uniqueN(x$sample), '\n')
