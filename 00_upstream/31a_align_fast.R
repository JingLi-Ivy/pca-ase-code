#!/usr/bin/env Rscript
# Vectorised version of 08a: align CPGEA VCF REF/ALT to 1000G; swap GT 0<->1 where alleles are swapped; drop the rest.
# (columns of the 1000G table are named rpos/rref/ralt: a column called "r" shadows the table inside data.table [ ])
args <- commandArgs(TRUE); suppressPackageStartupMessages(library(data.table))
hdr <- system(sprintf("zcat %s | grep '^#'", args[1]), intern = TRUE)
v <- fread(cmd = sprintf("zcat %s | grep -v '^#'", args[1]), header = FALSE, sep = "\t", colClasses = "character")
rf <- fread(args[2], header = FALSE, col.names = c("rpos", "rref", "ralt"), colClasses = "character")
dup <- rf$rpos[duplicated(rf$rpos)]; rf <- rf[!(rpos %chin% dup)]
m <- match(v$V2, rf$rpos); ok <- !is.na(m) & nchar(v$V4) == 1 & nchar(v$V5) == 1
keep <- ok & v$V4 == rf$rref[m] & v$V5 == rf$ralt[m]; swap <- ok & v$V4 == rf$ralt[m] & v$V5 == rf$rref[m]
keep[is.na(keep)] <- FALSE; swap[is.na(swap)] <- FALSE
for (j in 10:ncol(v)) set(v, which(swap), j, chartr("01", "10", v[[j]][swap]))
v[swap, `:=`(V4 = rf$rref[m[swap]], V5 = rf$ralt[m[swap]])]; v <- v[keep | swap]; v[, V8 := "."]
writeLines(hdr, args[3]); fwrite(v, args[3], append = TRUE, sep = "\t", col.names = FALSE, quote = FALSE)
cat(sprintf("keep %d swap %d drop %d\n", sum(keep), sum(swap), sum(!(keep | swap))))
