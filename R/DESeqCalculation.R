#' Differential Expression Analysis Using DESeq2
#'
#' @description
#' A wrapper function for DESeq2 analysis tailored for RNA-seq exon count data. 
#' The function executes the DESeq2 workflow, including size factor estimation, 
#' normalization, and differential expression testing between two conditions, such as
#' wild type (WT) and knockout (KO). It produces a series of diagnostic and 
#' exploratory plots to visualize data quality and expression differences, including 
#' boxplots, principal component analysis (PCA) plots, multi-dimensional scaling (MDS) 
#' plots, MA plots, histograms of p-values, and heatmaps for top differentially expressed 
#' genes.
#'
#' @param dat A matrix with RNAseq exon count data where rows are genes and columns 
#' are samples.
#' @param genotypes A factor vector indicating the genotype of each sample, used 
#' in the design formula for DESeq2.
#' @param fc The fold change threshold for determining significant differential 
#' expression. Default is 1.15.
#'
#' @return A list containing the DESeq2 results object, various ggplot and base 
#' R plot objects for visualization, and a pheatmap object for heatmap visualization.
#' @export
#' @examples
#' \dontrun{
#' # Simulate RNA-seq count data
#' set.seed(123) # for reproducibility
#' counts <- matrix(rpois(2000, lambda = 10), ncol = 10)
#' colnames(counts) <- paste0("Sample", 1:10)
#' rownames(counts) <- paste0("Gene", 1:200)
#' genotypes <- factor(rep(c("WT", "KO"), each = 5))
#'
#' # Perform the DESeq2 analysis
#' deseq_results <- DESeqCalculation(dat = counts, genotypes = genotypes)
#' 
#' # Explore the results
#' print(deseq_results$plots$PCAplot)
#' }
#' # Generate arbitrary toy data
#' dat <- matrix(rnbinom(n=5000, mu=100, size= 1/0.5),ncol=10)
#' for (x in 1:250) {
#' dat[x,] <- dat[x,] + sample(800:1200, 10, replace=TRUE)
#' }
#' for (x in 250:500) {
#' dat[x,] <- dat[x,] + sample(0:100, 10, replace=TRUE)
#' }
#' for (x in 1:5) {
#' dat[,x] <- dat[,x] + sample(800:1000, 5, replace=TRUE)
#' }
#' dat <- matrix(as.numeric(dat), ncol = ncol(dat))
#' colnames(dat) <- c("MeCP2_WT_1",	"MeCP2_WT_2",	"MeCP2_WT_3",
#' "MeCP2_WT_4",	"MeCP2_WT_5",	"MeCP2_KO_1",	"MeCP2_KO_2",	"MeCP2_KO_3",
#' "MeCP2_KO_4",	"MeCP2_KO_5")
#' vec <- as.character(seq(1,500, by = 1))
#' row.names(dat) <- vec
#' genotypes <- factor(c(rep("WT", 5), rep("KO", 5)), levels = c("KO", "WT"))
#'
#' # Run Differential Analysis
#' DESeqCalculation(dat = dat,genotypes = genotypes, fc = 1.15)
DESeqCalculation <- function(dat, genotypes, fc = 1.15){
  colData <- data.frame(samples = colnames(dat), genotypes = factor(genotypes),
                        row.names = colnames(dat))
  dds <- DESeqDataSetFromMatrix(countData = dat, colData = colData,
                                design = ~ genotypes)
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized=TRUE)
  print(dim(dat))

  ## size factor vs read counts
  sizeFactorRatio <- sizeFactors(dds)
  readCounts <- colSums(counts(dds))
  dds.mat <- data.frame(x = sizeFactors(dds), y = colSums(counts(dds))/1e6)
  p1 <- ggplot(dds.mat, aes(x = x, y = y)) +
    geom_point(aes(colour = dds$genotypes), size = 5) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey30") +
    xlab("Size Factor") +
    ylab("Number of Aligned Reads (in million)") + theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold", color = "black"),
          axis.text.y = element_text(size = 16, face = "bold", color = "black"),
          plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"),
          legend.position = c(.9, .1),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face = "bold", color = "black"))

  ## plots
  boxPlot(data = log2(dat+1), samples = factor(genotypes))
  rld.dds <- assay(rlog(dds, blind=FALSE))
  p2 <- PCAplot(data = rld.dds, genotypes = colData$genotypes,
                conditions = colData$genotypes)
  p2 <- p2 + ggrepel::geom_text_repel(aes(label = colData$samples),
                      size = 6, fontface = "bold",
                      color = "black", box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.2, "lines")) + theme_bw() +
    theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  p3 <- PCAplot(data = log2(dat+1), genotypes = colData$genotypes,
                conditions = colData$genotypes)
  p3 <- p3 + ggrepel::geom_text_repel(aes(label = colData$samples),
                      size = 6, fontface = "bold",
                      color = "black", box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.2, "lines")) + theme_bw() +
    theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  p4 <- MDSplot(data = log2(dat+1), genotypes = colData$genotypes,
                conditions = colData$genotypes)
  p4 <- p4 + ggrepel::geom_text_repel(aes(label = colData$samples), size = 6,
                      fontface = "bold",
                      color = "black", box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.2, "lines")) + theme_bw() +
    theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

  ## Expression Test
  dds <- DESeq(dds, betaPrior = TRUE)
  if(sum(unique(genotypes) %in% "KO") > 0){
    res.dds <- results(dds,contrast = c("genotypes","KO","WT"))
  }
  else{
    res.dds <- results(dds,contrast = c("genotypes","MUT","WT"))
  }
  res.dds$norm.counts <- counts(dds, normalized=TRUE)
  message(sum(res.dds$padj < 0.05, na.rm = TRUE))
  results <- as.data.frame(res.dds)

  ## Sorted as per adj. P-value
  resSort <- res.dds[order(res.dds$padj),]
  topGenes <- resSort[1:20,]
  topGenes$genes <- rownames(topGenes)

  ## Histogram and MA Plot with top 10 genes
  par(mfrow=c(1,1))
  hist(res.dds$pvalue[res.dds$baseMean > 1], breaks=0:20/20, col="grey50",
       border="white",
       main="Histogram of p-values with baseMean > 1")
  DESeq2::plotMA(resSort, main = "MA Plot")
  for(i in 1:nrow(topGenes)){
    with(topGenes[i, ],{
      points(topGenes$baseMean, topGenes$log2FoldChange, col="dodgerblue",
             cex=2, lwd=2)
      text(topGenes$baseMean, topGenes$log2FoldChange, topGenes$genes,
           pos=2, col="dodgerblue")})}

  ## Volcano Plot
  results <- results %>% tibble::rownames_to_column()
  colnames(results)[1] <- "gene"
  results <- results[which(!is.na(results$padj)),]
  results$Significant <- ifelse(results$log2FoldChange > log2(fc) &
                                  results$padj < 0.05, "Up",
                                ifelse(results$log2FoldChange < log2(1/fc)
                                & results$padj < 0.05, "Down","Not Signif"))
  pval.sig <- max(results[which(results$padj < 0.05),"pvalue"])
  p5 <- ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = c("green", "grey", "red")) +
    theme_bw(base_size = 16) +
    xlab("Log2 Fold Change") + ylab("-Log10 P-value") +
    geom_hline(aes(yintercept = -log10(pval.sig)), color="dodgerblue",
               linetype="dashed") +
    geom_vline(aes(xintercept = log2(fc)), color="dodgerblue",
               linetype="dashed") +
    geom_vline(aes(xintercept = log2(1/fc)), color="dodgerblue",
               linetype="dashed") +
    theme(axis.title = element_text(size = 22, face = "bold"),
        axis.text.x = element_text(size = 22, face = "bold", color = "black"),
        axis.text.y = element_text(size = 22, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

  ## normalized values for the heatmaps
  topVarGenes <- order(rowVars(rld.dds),decreasing=TRUE)
  mat <- rld.dds[topVarGenes[1:50], ]
  topVarGenes <- order(rowVars(rld.dds),decreasing=TRUE)
  annot <- data.frame(genotypes = factor(genotypes), row.names = colnames(dat))
  if(sum(levels(annot$genotypes) %in% c("WT", "MUT", "KO")) >= 2){
    annot$genotypes <- relevel(annot$genotypes, ref = "WT")
  }
  p1.top <- pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                     show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 9, legend = TRUE, filename = NA,
                     fontsize_col = 10, scale = "row", #fontface="bold",
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "euclidean",
                     annotation_col = annot)

  ## upregulated and downregulated genes
  ind.up <- which(results$log2FoldChange > log2(fc) & results$padj < 0.05)
  ind.down <- which(results$log2FoldChange < log2(1/fc) & results$padj < 0.05)
  up.reg <- results[ind.up, ]
  down.reg <- results[ind.down, ]

  ## upregulated heatmaps
  up.rld.dds <- rld.dds[rownames(rld.dds) %in% up.reg$gene,]
  up.topVarGenes <- order(rowVars(up.rld.dds),decreasing=TRUE)
  up.mat <- up.rld.dds[up.topVarGenes[1:50], ]
  p2.up <- pheatmap(up.mat, cluster_rows = TRUE, cluster_cols = TRUE,
                    show_rownames = TRUE, show_colnames = TRUE,
                    fontsize_row = 9, legend = TRUE, filename = NA,
                    fontsize_col = 10, scale = "row", #fontface="bold",
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "euclidean",
                    annotation_col = annot)

  ## downregulated heatmaps
  down.rld.dds <- rld.dds[rownames(rld.dds) %in% down.reg$gene,]
  down.topVarGenes <- order(rowVars(down.rld.dds),decreasing=TRUE)
  down.mat <- down.rld.dds[down.topVarGenes[1:50], ]
  p3.down <- pheatmap(down.mat, cluster_rows = TRUE, cluster_cols = TRUE,
                      show_rownames = TRUE, show_colnames = TRUE,
                      fontsize_row = 9, legend = TRUE, filename = NA,
                      fontsize_col = 10, scale = "row", #fontface="bold",
                      clustering_distance_rows = "correlation",
                      clustering_distance_cols = "euclidean",
                      annotation_col = annot)

  ## PCA plot on high variable genes
  p6 <- PCAplot(data = log2(dat[rowMeans(dat) > 30, ]+1),
                genotypes = colData$genotypes)

  print(p6)

  ## DEGs plot
  results.list <- list(results = results, up.reg = up.reg, down.reg = down.reg,
                      counts = dat, plots = list(sizefactorPlot = p1,
                      PCAplot_rld = p2, PCAplot = p3, MDSplot = p4,
                      Volplot = p5, PCAplot_highvar = p6, heatmap1 = p1.top,
                      heatmap2 = p2.up, heatmap3 = p3.down))
  return(results.list)
}
