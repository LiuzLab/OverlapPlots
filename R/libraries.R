#' @import  DESeq2  edgeR ggplot2 ggrepel gridExtra patchwork rafalib
#' @import reshape2 scales tibble pheatmap graphics matrixStats stats
# SummarizedExperiment cowplot dplyr
#' @importFrom grDevices palette
#' @importFrom utils globalVariables
#' @importFrom graphics points text
#' @importFrom gridExtra combine
#' @importFrom matrixStats rowRanges count
#' @importFrom patchwork align_plots
#' @importFrom stats start filter lag end
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr inner_join
#' @importFrom cowplot plot_grid
#' @importFrom details details
#' @importFrom sessioninfo session_info
#' @importFrom OpenRepGrid randomWords
#' @importFrom utils read.table
rafalib::mypar()
globalVariables(c("x","y","log2FoldChange", "Significant",
                  "samples", "sampleName", "normCounts",
                  "genotype", "PC1", "PC2", "pval.log10", "logFC.crude",
                  "gene.length", "refseq", "gene.type", "variable",
                  "value", "bin", "logFC"))






