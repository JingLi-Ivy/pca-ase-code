# Align REF/ALT of the CPGEA VCF to 1000G: swap GT 0<->1 when alleles are swapped; drop inconsistent or non-biallelic SNPs
args <- commandArgs(TRUE); suppressPackageStartupMessages(library(data.table))
con <- gzfile(args[1]); L <- readLines(con); close(con)
hdr <- L[startsWith(L, "#")]; body <- L[!startsWith(L, "#")]
ref <- fread(args[2], header = FALSE, col.names = c("pos","r","a"))
out <- character(0); n <- c(keep = 0, swap = 0, drop = 0)
for (ln in body) {
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]; pos <- as.integer(f[2]); p0 <- pos; r1 <- ref[ref$pos == p0]
  if (nrow(r1) != 1 || nchar(f[4]) != 1 || nchar(f[5]) != 1) { n["drop"] <- n["drop"] + 1; next }
  if (f[4] == r1$r && f[5] == r1$a) { n["keep"] <- n["keep"] + 1 }
  else if (f[4] == r1$a && f[5] == r1$r) { n["swap"] <- n["swap"] + 1; f[4] <- r1$r; f[5] <- r1$a
    g <- f[10:length(f)]; g <- chartr("01", "10", g); f[10:length(f)] <- g }
  else { n["drop"] <- n["drop"] + 1; next }
  f[8] <- "."; out <- c(out, paste(f, collapse = "\t")) }
writeLines(c(hdr, out), args[3]); print(n)
