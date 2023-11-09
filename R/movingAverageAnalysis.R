#' Comprehensive Analysis of RNA-seq Data with Moving Average Visualization
#'
#' @description
#' Integrates several steps of RNA-seq data analysis for knockout (KO) versus wild type (WT) comparisons
#' using the DESeq2 package. The function processes exon count data, executes differential expression 
#' analysis, performs k-means clustering, and generates overlap plots showing log2 fold changes 
#' versus gene length. Additionally, scatter plots for gene distribution and other relevant plots 
#' are produced.
#'
#' @param dat Path to the whole cell RNAseq exon counts file in .txt format.
#' @param genotypes A factor indicating the classification of samples, typically "WT" versus "KO".
#' @param degs.dat Path to the whole cell RNAseq differentially expressed genes (DEGs) file in .txt format.
#' @param bin.size Integer specifying the bin size for averaging in the overlap plots.
#' @param shift.size Integer specifying the step size for the moving average calculation in the overlap plots.
#' @param title A character string to be used as the title for the plots, indicating the dataset.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{res}{Object containing the results from DESeq2 analysis.}
#'   \item{overlapPlots}{List of ggplot objects for overlap plots, including combined and individual plots for all genes, edgeR DEGs, and DESeq2 DEGs.}
#'   \item{scatterPlot1}{Scatter plot for distribution of genes with log2 fold change greater than 0.}
#'   \item{scatterPlot2}{Scatter plot for distribution of genes with log2 fold change greater than log2(1.2).}
#' }
#' @export
#' @examples
#' \dontrun{
#' library(OverlapPlots)
#' data(countsfile)
#' data(refseq)
#' data(degsfile)
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' 
#' # Assuming 'countsfile', 'refseq', and 'degsfile' are loaded datasets
#' analysis_results <- movingAverageAnalysis(dat = countsfile, genotypes = genotypes, 
#'                                          degs.dat = degsfile, bin.size = 60, 
#'                                          shift.size = 6, 
#'                                          title = "KO/WT whole cell dataset")
#' print(analysis_results$overlapPlots$combined)
#' }
movingAverageAnalysis <- function(dat, genotypes,  degs.dat, bin.size, shift.size,
                     title = "MeCP2 KO"){
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

#' Methylation Context Analysis with Differential Gene Expression
#'
#' @description
#' Performs an integrated analysis of methylation context (mCA) data with RNA-seq 
#' differential expression results. It generates comparative plots to examine log2 
#' fold change versus methylation context (mC), distribution of gene lengths within 
#' bins, and the contribution of methylation from long versus short genes.
#'
#' @param count.file A dataframe resulting from the `movingAverageAnalysis` containing 
#' normalized count data along with gene names and fold changes.
#' @param degs.dat Path to the file with differential gene expression data from RNA-seq 
#' analysis in .txt format.
#' @param mCA A dataframe with methylation context data, including gene names, mC values, 
#' and gene lengths.
#' @param bin.size Integer specifying the bin size for averaging in the plots.
#' @param shift.size Integer specifying the step size for the moving average calculation 
#' in the plots.
#' @param title A string indicating the title of the dataset for plot annotations. 
#' Default is "MeCP2 KO".
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{res}{An object with various analysis results and data subsets used in plotting.}
#'   \item{combined}{A composite ggplot object that includes all the generated plots.}
#' }
#' @export
#' @examples
#' \dontrun{
#' library(OverlapPlots)
#' # Assuming 'countsfile', 'degsfile', and 'mCA' are loaded datasets
#' dat <- countsfile # previously loaded or created dataset
#' degs.dat <- degsfile # previously loaded or created DEGs dataset
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' mCA.sub <- mCA[mCA$CA >= 5 & mCA$gene.length >= 4500,]
#' mCA.sub$mCA.CA <- mCA.sub$mCA/mCA.sub$CA
#'
#' # Perform moving average analysis
#' wholeCell.KO <- movingAverageAnalysis(dat, genotypes, degs.dat, bin.size = 60, 
#'                                       shift.size = 6, title = "KO/WT whole cell dataset")
#' # Perform mCA analysis
#' mCA.results <- movingAverageAnalysismCA(wholeCell.KO$res$results, degs.dat, mCA.sub,
#'                                         bin.size = 60, shift.size = 6, 
#'                                         title = "KO/WT whole cell dataset")
#' print(mCA.results$combined)
#' }


movingAverageAnalysismCA <- function(count.file, degs.dat, mCA, bin.size, shift.size,
                        title = "MeCP2 KO"){

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
