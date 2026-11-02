suppressPackageStartupMessages(library(data.table))
fs <- Sys.glob("/path/to/data_root/ASE/ADASTRA6/ADASTRA.v.6.1.Mabel/TF/*_HUMAN.tsv")
mp <- fread("/path/to/data_root/ASE/code/adastra_tf_map.csv")
out <- rbindlist(lapply(fs, function(f) { l <- system(sprintf("grep -P \"^chr17\\t\\d+\\t44947872\\t\" %s", f), intern = TRUE)
  if (!length(l)) return(NULL); h <- strsplit(readLines(f, 1), "\t")[[1]]; h[1] <- "chr"
  d <- fread(text = paste(l, collapse = "\n"), header = FALSE, sep = "\t", col.names = h, colClasses = "character")
  d[, mnem := sub("_HUMAN.tsv", "", basename(f))]; d }), fill = TRUE)
num <- c("es_mean_ref","es_mean_alt","fdrp_bh_ref","fdrp_bh_alt","motif_fc","n_aggregated","total_cover"); out[, (num) := lapply(.SD, as.numeric), .SDcols = num]
print(names(mp)); out <- merge(out, mp, by.x = "mnem", by.y = "mn", all.x = TRUE)
out[, tf := fifelse(is.na(sym) | sym == "", mnem, sym)]
out[, `:=`(pref = fifelse(fdrp_bh_ref < fdrp_bh_alt, "G", "C"), minfdr = pmin(fdrp_bh_ref, fdrp_bh_alt))]
fwrite(out[, .(tf, mnem, n_aggregated, total_cover, es_mean_ref, es_mean_alt, fdrp_bh_ref, fdrp_bh_alt, motif_fc, pref, minfdr)][order(minfdr)], "adastra_rs9914422.tsv", sep = "\t")
print(nrow(out)); print(out[order(minfdr)][1:12, .(tf, mnem, es_mean_ref, es_mean_alt, pref, minfdr, motif_fc)]); print(out[tf %in% c("FOXM1","AR","FOXA1","MYBL2","E2F4")])
