
count.file <- "dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt"
countsfile <- read.table(file = count.file, sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
usethis::use_data(countsfile, overwrite = FALSE)
