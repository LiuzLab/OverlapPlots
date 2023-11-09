#' Overlap Plot Wrapper for Gene Expression Analysis
#'
#' @description
#' This function acts as a wrapper for generating overlap plots from differential
#' expression analysis results. It enriches the dataset with reference sequence
#' information, computes log2 fold changes between different groups, and creates
#' running average plots (log2 fold change vs mean gene length). It utilizes the
#' `overlayGabelsPlot` function for plotting, which internally calculates running
#' averages and p-values from a student's t-test, and then produces overlay line
#' plots with confidence intervals.
#'
#' @param dat A matrix or dataframe containing the differential expression analysis
#' results from `DESeqCalculation`. It should include gene names and log fold change
#' values among other analysis results.
#' @param refseq A dataframe containing reference sequence information that will be
#' joined with `dat` on the gene names.
#' @param KO.idx Indices or column names in `dat` corresponding to the treatment group.
#' @param WT.idx Indices or column names in `dat` corresponding to the control group.
#' @param WT1.idx Subset indices within the control group for the first comparison.
#' @param WT2.idx Subset indices within the control group for the second comparison.
#' @param bin.size The size of the window (number of genes) over which to calculate
#' the running average.
#' @param shift.size The number of genes to shift the window by for each calculation
#' of the running average.
#' @param confidenceinterval The width of the confidence interval to be displayed on
#' the plots. Defaults to 0.50 (50% confidence interval).
#' @param shrink_lfc Logical indicating whether to use shrunken log2 fold changes. 
#' Defaults to `FALSE`.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{res}{A list returned from `overlayGabelsPlot` containing ggplot objects for the overlay plots.}
#'   \item{plot}{A combined ggplot object with the running average plot and the p-value plot.}
#'   \item{log2FC.length}{A dataframe with the gene names, log2 fold change, and gene length used for the plot.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' load(file = "path/to/mm10_ncbi-refSeqGene_Dec2019.RData")
#' genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
#' dat <- as.data.frame(matrix(runif(1000), ncol = 10))
#' colnames(dat) <- c("gene1", "gene2", ..., "gene10")
#' # Assume DESeqCalculation returns a dataframe like structure
#' wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
#' mat <- wholeCell.KO$results[, c(8:27, 1, 3)]
#' colnames(mat)[21] <- "gene.name"
#' grp.idx <- WTgrpKmeans(control_mat = mat[, 1:10])
#' res <- overlapWrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
#'                       WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
#'                       WT2.idx = grp.idx$WT.idx2, bin.size = 200,
#'                       shift.size = 40, confidenceinterval = 0.50)
#' print(res$plot)
#' }
overlapWrapper <- function(dat, refseq, KO.idx, WT.idx, WT1.idx, WT2.idx,
                            bin.size, shift.size, confidenceinterval = 0.50,
                           shrink_lfc = FALSE){
  dat.annot <- inner_join(x = dat, y = refseq, by = "gene.name")
  cat("Control Group 1:",WT1.idx,"\n")
  cat("Control Group 2:",WT2.idx,"\n")
  gene.idx <- which(colnames(dat.annot) %in% "gene.name")
  log2FC.WT <- dat.annot[,c(WT.idx, gene.idx)]
  log2FC.WT$comp.mat <- apply(log2FC.WT[, WT.idx], 1,
                              function(r){log2((mean(r[WT1.idx]) + 1) /
                                                 (mean(r[WT2.idx]) +1))})
  if(shrink_lfc == FALSE){
    log2FC.KO <- logofMeansBetweenAB(dat = dat.annot, A.samples = WT.idx,
                                        B.samples = KO.idx)
    log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
                                y = log2FC.KO[,c("gene.name","logFC.crude",
                                                 "gene.length")],
                                by = "gene.name")
  }else{
    log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
                                y = dat.annot[,c("gene.name",
                                                 "log2FoldChange",
                                                 "gene.length")],
                                by = "gene.name")
  }

  res <- overlayGabelsPlot(mat = log2FC.length[,c(2:4)],
                             comp.between1 = "(WT/WT)",
                             comp.between2 = "(KO/WT)", bin.size, shift.size,
                           confidenceinterval)
  plot.margin <- unit(c(1,0.5,0.5,0.5), "cm")
  res.plot <- plot_grid(res$plot1 + coord_cartesian(ylim = c(-0.4,0.4)),
                        res$plot2 + coord_cartesian(ylim = c(0,50)), ncol = 1,
                        align = 'v') + theme(plot.margin = plot.margin)
  return(list(res = res, plot = res.plot, log2FC.length = log2FC.length))
}

#' Perform K-Means Clustering on Control Group Data
#'
#' @description
#' Performs k-means clustering on a control (WT) subset from differential expression analysis results.
#' This function is designed to identify patterns or groups within control samples based on their gene
#' expression profiles. It uses the standard k-means clustering algorithm to partition samples into a
#' specified number of clusters and returns the indices of samples in each cluster.
#'
#' @param control_mat A matrix or dataframe containing the control (WT) subset from the
#' differential expression analysis results. Rows correspond to genes and columns to samples.
#' @param centers The desired number of clusters (centers) to find within the control data.
#' Defaults to 2.
#' @param iter.max The maximum number of iterations allowed for the k-means algorithm.
#' Defaults to 1000.
#'
#' @return A list with two elements, each containing the indices of the control samples
#' assigned to cluster 1 and cluster 2, respectively. The list elements are named `WT.idx1`
#' and `WT.idx2`.
#' @export
#'
#' @examples
#' # Simulate gene expression data for a control group with two distinct patterns
#' set.seed(123) # For reproducibility
#' control_data <- rbind(
#'   matrix(rnorm(500, sd = 0.3), ncol = 10),
#'   matrix(rnorm(500, mean = 1, sd = 0.3), ncol = 10)
#' )
#' colnames(control_data) <- paste0("Sample", 1:10)
#' rownames(control_data) <- paste0("Gene", 1:100)
#' 
#' # Perform k-means clustering
#' kmeans_results <- WTgrpKmeans(control_mat = control_data)
#' 
#' # kmeans_results$WT.idx1 contains indices of samples in cluster 1
#' # kmeans_results$WT.idx2 contains indices of samples in cluster 2

WTgrpKmeans <- function(control_mat, centers = 2, iter.max = 1000){
  group <- kmeans(t(control_mat), centers, iter.max)$cluster
  idx1 <- which(group %in% 1)
  idx2 <- which(group %in% 2)
  return(list(WT.idx1 = idx1, WT.idx2 = idx2))
}

#' Equal Size K-Means Clustering on Control Group Data
#'
#' @description
#' Performs a variation of k-means clustering on the control (WT) subset from 
#' differential expression analysis results to create clusters of approximately 
#' equal size. This method modifies the standard k-means approach by ordering 
#' points based on their distance to the closest center minus the distance to 
#' the farthest cluster and then assigns each point to the best cluster in this order.
#' The method aims to address the challenge of unequal cluster sizes often 
#' encountered in standard k-means.
#'
#' @param control_mat A matrix or dataframe containing the control (WT) subset 
#' from differential expression analysis results. Rows correspond to genes and 
#' columns to samples.
#' @param centers The desired number of clusters to find within the control data. 
#' Defaults to 2.
#' @param iter.max The maximum number of iterations allowed for the k-means 
#' algorithm. Defaults to 1000.
#'
#' @return A list containing two elements: `WT.idx1` and `WT.idx2`, which are 
#' indices of the control samples assigned to clusters 1 and 2, respectively.
#' @export
#'
#' @examples
#' # Simulate gene expression data for a control group
#' set.seed(123) # For reproducibility
#' control_data <- rbind(
#'   matrix(rnorm(1000, sd = 0.3), ncol = 10),
#'   matrix(rnorm(1000, mean = 1, sd = 0.3), ncol = 10)
#' )
#' colnames(control_data) <- paste0("Sample", 1:10)
#' rownames(control_data) <- paste0("Gene", 1:100)
#'
#' # Perform equal size k-means clustering
#' kmeans_results <- WTgrpKmeansEqualSize(control_mat = control_data)
#' 
#' # kmeans_results$WT.idx1 contains indices of samples in cluster 1
#' # kmeans_results$WT.idx2 contains indices of samples in cluster 2

WTgrpKmeansEqualSize <- function(control_mat, centers = 2, iter.max = 1000){
  size <- ceiling(nrow(t(control_mat))/centers)
  group <- kmeans(t(control_mat), centers, iter.max)
  new_group <- rep(NA, nrow(t(control_mat)))
  new_centers <- lapply(1:centers, function(r){
    euc_dist <- control_mat - group$centers[r,]
    sqrt(apply(euc_dist, 2, function(x) sum(x^2)))})
  new_centers <- matrix(unlist(new_centers), ncol = centers)
  new_clust_size <- rep(0, centers)
  sample_ord <- order(apply(new_centers, 1, min) - apply(new_centers, 1, max))
  for(i in sample_ord){
    bestcl <- which.max(new_centers[i,])
    new_group[i] <- bestcl
    new_clust_size[bestcl] <- new_clust_size[bestcl] + 1
    if(new_clust_size[bestcl] >= size){
      new_centers[,bestcl] <- NA
    }
  }
  idx1 <- which(new_group %in% 1)
  idx2 <- which(new_group %in% 2)
  return(list(WT.idx1 = idx1, WT.idx2 = idx2))
}


#' Generate Overlap Plots for Differential Gene Expression and Methylation Analysis
#'
#' @description
#' This function creates overlap plots for differentially expressed genes (DEGs) with 
#' methylation data. It performs an inner join between DEGs data and reference sequence 
#' information, computes log2 fold changes between control groups, and then calls the 
#' `overlaymC` function to generate the plots, which illustrate log2 fold change versus 
#' mean gene length, including methylation context (mCA/CA).
#'
#' @param degs.dat A dataframe containing DEGs results from differential expression 
#' analysis, such as from `DESeqCalculation`. The dataframe should include gene names, 
#' log fold changes, gene lengths, and FDR values.
#' @param count.dat A dataframe with counts data from which log2 fold changes will be 
#' calculated. It should have gene names and count data for each sample.
#' @param refseq A dataframe with reference sequence information including methylation data 
#' (mCA/CA), to be joined with `degs.dat`.
#' @param WT1.idx Indices of the first control group samples in `count.dat`.
#' @param WT2.idx Indices of the second control group samples in `count.dat`.
#' @param bin.size The size of the bin for calculating running averages in the plot.
#' @param shift.size The step size for the moving window when calculating running averages.
#' @param methyl.type A string indicating the type of methylation context, typically "mCA/CA".
#' @param degs Logical flag indicating whether to filter `degs.dat` by FDR < 0.05 
#' to consider only significant DEGs. Defaults to `TRUE`.
#'
#' @return A list containing the overlay plot object generated by `overlaymC` and other 
#' relevant data for further analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' load(file = "path/to/mm10_ncbi-refSeqGene_Dec2019.RData")
#' mCA <- data.frame(readRDS("path/to/CAperGene_10wk_boxer.RDS"), stringsAsFactors = FALSE)
#' mCA <- mCA[mCA$CA >= 5 & mCA$gene.length >= 4500,]
#' mCA$mCA.CA <- mCA$mCA/mCA$CA
#' mCA1 <- mCA[, c(1, 6, 7)]
#' 
#' degs.dat <- read.table("path/to/KO-WT_whole-cell_RNA-seq.txt", 
#'                        sep = "\t", header = TRUE, row.names = 1, 
#'                        stringsAsFactors = FALSE)
#' mat <- wholeCell.KO$results[, c(8:27, 1, 3)]
#' colnames(mat)[21] <- "gene.name"
#' 
#' grp.idx <- WTgrpKmeans(control_mat = mat[, 1:10])
#' if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
#'   message("Cluster size is not equal, running equal size k-means variation.")
#'   grp.idx <- WTgrpKmeansEqualSize(control_mat = mat[, 1:10])
#' }
#' 
#' res <- overlapDegsmCAWrapper(degs.dat = degs.dat, count.dat = mat, 
#'                              refseq = mCA1, WT1.idx = grp.idx$WT.idx1, 
#'                              WT2.idx = grp.idx$WT.idx2, bin.size = 60, 
#'                              shift.size = 6, methyl.type = "mCA/CA")
#' print(res$plot)
#' }

overlapDegsmCAWrapper <- function(degs.dat, count.dat, refseq, WT1.idx,
                                     WT2.idx, bin.size, shift.size, methyl.type,
                                     degs = TRUE){
  if(degs == TRUE){
    degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
  }
  cat("Number of genes =",dim(degs.dat)[1],"\n\n")
  degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
                         refseq, #%>% rownames_to_column(var = "gene.name"),
                         by = "gene.name")
  ## overlap plot
  cat("Control Group 1:",WT1.idx,"\n")
  cat("Control Group 2:",WT2.idx,"\n")
  log2FC.WT <- data.frame(count.dat[,c(1:10, 21)], stringsAsFactors = FALSE)
  log2FC.WT$comp.mat <- apply(count.dat[,c(WT1.idx, WT2.idx)], 1,
                              function(r){log2((mean(r[1:length(WT1.idx)])+1) /
                                                 (mean(r[(length(WT1.idx)+1):10])+1))})
  log2FC.length <- inner_join(x = log2FC.WT[,c("gene.name","comp.mat")],
                              y = degs.dat[,c("gene.name","logFC",
                                              "gene.length","mCA.CA")],
                              by = "gene.name")
  message(dim(log2FC.length)[1])
  print(log2FC.length[,c(2:5,1)])
  res <- overlaymC(mat = log2FC.length[,c(2:5,1)], comp.between1 = "(WT/WT)",
                    comp.between2 = "(KO/WT)", bin.size = bin.size,
                    shift.size = shift.size, methyl.type = methyl.type)
  return(res = res)
}

