#' Global Main analysis function
#'
#' @description
#' This function is intended to generate a figure where it includes the Overlap
#' plots comparing log2FC vs the gene length. This plot is generated using Boxer
#'  et al KO/WT whole cell dataset.
#' @param dat whole cell RNAseq exon counts counts .txt format
#' @param genotypes factor WT vs KO
#' @param degs.dat whole cell RNA seq DEGs .txt file
#' @param bin.size bin size integer
#' @param shift.size shift size integer
#' @param title dataset title
#'
#' @return main analysis plots
#' @export
#' @examples
#' library(OverlapPlots)
#' data(countsfile)
#' data(refseq)
#' data(degsfile)
#' dat <- countsfile
#' degs.dat <- degsfile
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' bin.size <- 60
#' shift.size <- 6
#' wholeCell.KO <- movingAverageAnalysis(dat, genotypes,  degs.dat, bin.size, shift.size,
#'                         title = "KO/WT whole cell dataset")

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

#' Global Main analysis function for mCA
#'
#' @description
#' This figure illustrates the comparison between the Log2FC vs the mC (Fig A).
#' It further looks into the distribution of the long and short genes in each
#' bin (Fig B and Fig C) and then investigate if contribution of higher mC comes
#' from long genes as compared to the short genes (Fig D).
#'
#' @param count.file movingAverageAnalysis results
#' @param degs.dat whole cell RNA seq DEGs .txt file
#' @param mCA mCA data
#' @param bin.size bin size integer
#' @param shift.size shift size integer
#' @param title dataset title
#'
#' @return main analysis plots for mCA
#' @export
#' @examples
#' library(OverlapPlots)
#' data(countsfile)
#' data(refseq)
#' data(degsfile)
#' data(mCA)
#' dat <- countsfile
#' degs.dat <- degsfile
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' bin.size <- 60
#' shift.size <- 6
#' wholeCell.KO <- movingAverageAnalysis(dat, genotypes,  degs.dat, bin.size, shift.size,
#'                         title = "KO/WT whole cell dataset")
#' mCA <- mCA[mCA$CA >= 5,]
#' mCA <- mCA[mCA$gene.length >= 4500,]
#' mCA$mCA.CA <- mCA$mCA/mCA$CA
#' mCA.sub <- mCA[,c(1,6,7)]
#' mCA.wholeCell.KO <- movingAverageAnalysismCA(wholeCell.KO$res$results, degs.dat, mCA.sub,
#'                                bin.size, shift.size,
#'                                title = "KO/WT whole cell dataset")

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
