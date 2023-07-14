# # Creating toy examples for exported functions
#
# # Log means
# dat <- matrix(rnorm(600, sd = 0.3), ncol = 6)
# dat <- as.data.frame(dat)
# logofMeans.between.A.B(dat = dat, A.samples = 1:3, B.samples = 4:6)
#
# # Gabel function
#
# # mat is dat.annotated (data counts + ref sequence) with log means results (log fold change and others...)
#
# gabels.plot <- function(mat, length.type = "Gene", comp.between = "",
#                         y.axis = "Mean Log2 Fold Change"){
#   p1 <- moving.average.function(dat = mat, bin.size = 200, shift.size = 40,
#                                 length.type, comp.between, y.axis)
#   plot(p1)
# }

# Wrapper Need refseq information toy data

## Load annotation dataset
load(file = "../../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))

## Load gene counts
dat <- read.table("../../dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)

# Run DESeq analysis
wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
mat <- wholeCell.KO$results[,c(8:27,1,3)]
colnames(mat)[21] <- "gene.name"

## same cluster sizes
grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])

## Generate overlap plots
res1 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 200,
                        shift.size = 40)
res1$plot

#' \dontrun{
#' # Generate arbitrary toy data
#' dat <- matrix(rnbinom(n=10000, mu=100, size= 1/0.5),ncol=20)
#' for (x in 1:250) {
#' dat[x,] <- dat[x,] + sample(800:1200, 10, replace=T)
#' }
#' for (x in 250:500) {
#' dat[x,] <- dat[x,] + sample(0:100, 10, replace=T)
#' }
#' for (x in 1:5) {
#' dat[,x] <- dat[,x] + sample(800:1000, 5, replace=T)
#' }
#' dat <- matrix(as.numeric(dat), ncol = ncol(dat))
#' colnames(dat) <- c("MeCP2_WT_1",	"MeCP2_WT_2",	"MeCP2_WT_3",
#' "MeCP2_WT_4",	"MeCP2_WT_5",	"MeCP2_WT_6", "MeCP2_WT_7", "MeCP2_WT_8",
#' "MeCP2_WT_9", "MeCP2_WT_10", "MeCP2_KO_1",	"MeCP2_KO_2",	"MeCP2_KO_3",
#' "MeCP2_KO_4",	"MeCP2_KO_5", "MeCP2_KO_6", "MeCP2_KO_7", "MeCP2_KO_8",
#' "MeCP2_KO_9", "MeCP2_KO_10")
#' vec <- as.character(seq(1,500, by = 1))
#' row.names(dat) <- vec
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' # Run Differential Analysis
#' wholeCell.KO <- DESeqCalculation(dat = dat,genotypes = genotypes, fc = 1.15)
#' # Plot figures
#' mat <- wholeCell.KO$results[,c(8:27,1,3)]
#' colnames(mat)[21] <- "gene.name"
#' grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
#' res1 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
#'                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#'                         WT2.idx = grp.idx$WT.idx2, bin.size = 200,
#'                         shift.size = 40)
#' res1$plot
#' }
