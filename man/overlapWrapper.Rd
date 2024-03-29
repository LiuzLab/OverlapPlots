% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlap_wrappers.R
\name{overlapWrapper}
\alias{overlapWrapper}
\title{Overlap Plot Wrapper for Gene Expression Analysis}
\usage{
overlapWrapper(
  dat,
  refseq,
  KO.idx,
  WT.idx,
  WT1.idx,
  WT2.idx,
  bin.size,
  shift.size,
  confidenceinterval = 0.5,
  shrink_lfc = FALSE
)
}
\arguments{
\item{dat}{Differential analysis results from "DESeqCalculation". It's expected to be a data frame or matrix with rows 
representing genes and columns representing different samples or experimental conditions.}

\item{refseq}{A reference sequence dataset, which should align with the genes represented in the 'dat' parameter.}

\item{KO.idx}{Treatment group index specifying the columns in 'dat' corresponding to the Knock-Out (KO) group.}

\item{WT.idx}{Control group index specifying the columns in 'dat' corresponding to the Wild-Type (WT) group.}

\item{WT1.idx}{WT1 index for subgroup analysis within the control group.}

\item{WT2.idx}{WT2 index for another subgroup within the control group.}

\item{bin.size}{An integer specifying the bin size for the moving average calculation.}

\item{shift.size}{An integer specifying the shift size for the moving average calculation.}

\item{confidenceinterval}{The confidence interval for the t-test calculations, defaulting to 0.50.}

\item{shrink_lfc}{A boolean indicating whether to shrink log2 fold changes, defaulting to FALSE.}
}
\value{
A list containing the results of the analysis, the generated plot, and the log2 fold change data.
}
\description{
This function is a specialized wrapper for generating overlap plots in the context of gene expression analysis. 
It integrates DESeqCalculation results with k-means clustering index groups to compute running average plots, 
illustrating the relationship between log2 fold change and mean gene length.
}
\details{
The function employs Overlay Gabel's moving average technique to calculate the mean values of log2 fold change 
and genomic lengths. It also computes p-values using the student t-test. The final output includes overlay line plots 
that visually represent these computations.
}
\examples{
\dontrun{
# Loading necessary datasets
load(file = "../../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))

# Reading and preparing gene count data
dat <- read.table(paste("../../dat/counts/GSE128178_10WT_",
                        "10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt"),
                 sep = "\t", stringsAsFactors=FALSE, header = TRUE, row.names = 1)
# Running DESeq analysis
wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
mat <- wholeCell.KO$results[,c(8:27,1,3)]
colnames(mat)[21] <- "gene.name"

# K-means clustering for subgroup analysis
grp.idx <- WTgrpKmeans(control_mat = mat[,1:10])

# Generating and visualizing overlap plots
res1 <- overlapWrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                       WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                       WT2.idx = grp.idx$WT.idx2, bin.size = 200,
                       shift.size = 40)
print(res1$plot)
}
}
