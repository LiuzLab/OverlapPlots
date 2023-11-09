#' Comprehensive Analysis Pipeline for Differential Expression Data
#'
#' @description
#' Conducts a complete analysis for RNA-seq data by integrating differential expression 
#' analysis (DESeq2), k-means clustering for grouping genes, and generating overlay plots 
#' to compare log2 fold changes versus gene length. The function processes count data, 
#' applies the DESeq2 workflow, performs clustering, reads in DEGs data, and creates a 
#' series of plots to visualize the findings.
#'
#' @param count.file Path to a file containing whole cell RNAseq exon counts in .txt format.
#' @param genotypes Factor indicating the classification of samples, typically "WT" versus "KO".
#' @param degs.file Path to a file containing whole cell RNAseq differentially expressed genes (DEGs) 
#' in .txt format.
#' @param bin.size Integer specifying the number of genes in each bin for calculating running averages.
#' @param shift.size Integer specifying the step size for the moving window when calculating running averages.
#' @param title A title for the dataset, used in plot annotations. Default is "MeCP2 KO".
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{res}{The DESeq2 analysis results object.}
#'   \item{overlapPlots}{A list of ggplot objects for overlay plots, including combined and individual plots 
#'                        for all genes, edgeR DEGs, and DESeq2 DEGs.}
#'   \item{scatterPlot1}{A ggplot object for the scatter plot with log2FC thresholded at 0.}
#'   \item{scatterPlot2}{A ggplot object for the scatter plot with log2FC thresholded at log2(1.2).}
#' }
#' @noRd
#'
#' @examples
#' \dontrun{
#' count.file <- "path/to/your/count_data.txt"
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' degs.file <- "path/to/your/degs_data.txt"
#' bin.size <- 60
#' shift.size <- 6
#' title <- "Sample Analysis Title"
#'
#' analysis_results <- finalRun(count.file = count.file, genotypes = genotypes, 
#'                              degs.file = degs.file, bin.size = bin.size, 
#'                              shift.size = shift.size, title = title)
#' print(analysis_results$overlapPlots$combined)
#' }
finalRun <- function(count.file, genotypes, degs.file, bin.size, shift.size,
                     title = "MeCP2 KO"){
  ## counts file
  dat <- read.table(file = count.file, sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE, row.names = 1)

  ## running DESeq and all the plots are in dds object
  dds <- DESeqCalculation(dat = dat, genotypes = genotypes, fc=1.15)
  mat <- dds$results[,c(8:27,1,3)]
  colnames(mat)[21] <- "gene.name"
  grp.idx <- WTgrpKmeans(control_mat = mat[,1:10])
  if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
    message("Cluster size is not equal, therefore run same size k-means
                variation!")
    grp.idx <- WTgrpKmeansEqualSize(control_mat = mat[,1:10])
  }

  ## overlap plot using all the genes
  op1 <- overlapWrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                         WT2.idx = grp.idx$WT.idx2, bin.size = 200,
                         shift.size = 40)

  ## Using the DEGs list
  degs.dat <- read.table(file = degs.file, sep = "\t", stringsAsFactors = FALSE,
                         header = TRUE, row.names = 1)
  degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
  cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
  degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
                         refseq, by = "gene.name")

  ## scatter plot (distribution between 2 different logFC)
  scater1 <- plotScatter(log2FC = log2(1), comp.between = "(> log2FC(0))",
                          dat = degs.dat[,c("gene.name","logFC", "FDR",
                                            "gene.length")])
  scater2 <- plotScatter(log2FC = log2(1.2), comp.between = "(> log2(1.2))",
                          dat = degs.dat[,c("gene.name","logFC", "FDR",
                                            "gene.length")],)

  ## prep for overlap plots with DEGs genes
  degs.dat <- inner_join(y = degs.dat[,c("gene.name", "logFC")],
                         x = dds$results[dds$results$gene %in%
                                           degs.dat$gene.name,c(8:27,1,3)],
                         by = c("gene" = "gene.name"))
  mat <- degs.dat[,c(1:21,23)]
  colnames(mat)[21] <- "gene.name"
  colnames(mat)[22] <- "log2FoldChange"

  ## logFC from degs from edgeR (Boxer et al.)
  op2 <- overlapWrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                         WT2.idx = grp.idx$WT.idx2, bin.size = bin.size,
                         shift.size = shift.size, shrink_lfc = TRUE)

  ## logFC from degs from DESeq2
  mat <- degs.dat[,c(1:22)]
  colnames(mat)[21] <- "gene.name"
  op3 <- overlapWrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                         WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                         WT2.idx = grp.idx$WT.idx2, bin.size = bin.size,
                         shift.size = shift.size, shrink_lfc = TRUE)
  ## plots
  p1 <- op1$plot + ggtitle("All genes") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  p2  <- op2$plot + ggtitle("edgeR DEGs") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  p3 <- op3$plot + ggtitle("DESeq2 DEGs") +
    theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
  q1 <- (p1 | p2 | p3) + plot_annotation(title = title, theme = theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")))

  ## return
  return(list(res = dds, overlapPlots = list(combined = q1, all_genes = op1,
                                             degs_edgeR = op2, degs_deseq2 = op3), scatterPlot1 = scater1,
              scatterPlot2 = scater2))
}


#' Main Analysis Function for Methylation Context (mC) Data
#'
#' @description
#' Performs a comprehensive analysis integrating methylation context with differential
#' gene expression data. The function generates multiple plots:
#' (A) Log2 fold change (Log2FC) versus methylation context (mC),
#' (B) Distribution of long and short genes within each bin,
#' (C) Comparative analysis of gene length distribution in bins,
#' (D) Examination of higher mC contribution from long genes versus short genes.
#'
#' @param count.file Path to a file containing whole cell RNAseq exon counts in 
#' .txt format. Expected to contain gene expression counts for different samples.
#' @param degs.file Path to a file containing differentially expressed genes (DEGs) 
#' from whole cell RNAseq analysis in .txt format.
#' @param mCA Dataframe containing methylation context data for genes.
#' @param bin.size Integer specifying the number of genes in each bin for calculating
#' running averages.
#' @param shift.size Integer specifying the step size for the moving window when 
#' calculating running averages.
#' @param title Title for the dataset to be used as an annotation in plots. 
#' Default is "MeCP2 KO".
#'
#' @return A list containing the results of the analysis and a combined ggplot object 
#' with the generated plots (A, B, C, D) arranged in a composite figure.
#' @noRd
#'
#' @examples
#' \dontrun{
#' count.file <- "path/to/your/whole_cell_counts.txt"
#' degs.file <- "path/to/your/whole_cell_degs.txt"
#' mCA_data <- readRDS("path/to/your/mca_data.RDS")
#' bin.size <- 60
#' shift.size <- 6
#' title <- "Sample Analysis Title"
#'
#' analysis_results <- finalRunmCA(count.file = count.file, degs.file = degs.file,
#'                                 mCA = mCA_data, bin.size = bin.size,
#'                                 shift.size = shift.size, title = title)
#' print(analysis_results$combined)
#' }
finalRunmCA <- function(count.file, degs.file, mCA, bin.size, shift.size,
                         title = "MeCP2 KO"){
  degs.dat <- read.table(file = degs.file, sep = "\t", stringsAsFactors = FALSE,
                         header = TRUE, row.names = 1)
  count.file <- count.file[,c(8:27,1,3)]
  colnames(count.file)[21] <- "gene.name"
  grp.idx <- WTgrpKmeans(control_mat = count.file[,1:10])
  if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
    message("Cluster size is not equal, therefore run same size k-means
                variation!")
    grp.idx <- WTgrpKmeansEqualSize(control_mat = count.file[,1:10])
  }
  res <- overlapDegsmCAWrapper(degs.dat = degs.dat, count.dat = count.file,
                                  refseq = mCA, WT1.idx = grp.idx$WT.idx1,
                                  WT2.idx = grp.idx$WT.idx2,
                                  bin.size = bin.size,
                                  shift.size = shift.size,
                                  methyl.type = "mCA/CA")
  comb.plot <- ((res$overlapPlot + coord_cartesian(ylim = c(-0.25,0.25))) /
                  res$mCAvGl_plot) | (res$plotBar / res$diffPlot)
  comb.plot <- comb.plot + plot_annotation(title = title, theme = theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")))
  return(list(res = res, combined = comb.plot))
}
