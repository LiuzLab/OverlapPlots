
degs.file <- "dat/DEGs/KO-WT_whole-cell_RNA-seq.txt"
degsfile <- read.table(file = degs.file, sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
usethis::use_data(degsfile, overwrite = TRUE)
