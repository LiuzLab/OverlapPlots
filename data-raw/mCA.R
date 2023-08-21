## code to prepare `mCA` dataset goes here

mCA <- data.frame(readRDS("dat/mCA/CAperGene_10wk_boxer.RDS"), stringsAsFactors = FALSE)

usethis::use_data(mCA, overwrite = TRUE)
