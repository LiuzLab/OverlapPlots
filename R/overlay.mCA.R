## overlay plot

#' Overlay Plots for Methylation and Gene Expression Analysis
#'
#' This function generates overlay plots to analyze the relationship between methylation levels 
#' (mCA/CA) and gene expression across different gene lengths. It uses binning to aggregate data points 
#' and visualize trends in methylation and expression changes.
#'
#' @param mat A dataframe with columns corresponding to comparison matrix values (comp.mat),
#' log fold change (logFC.crude), gene length, methylation values (mCA.CA), and gene names.
#' @param bin.size Numeric, the size of the bin for aggregating genes.
#' @param shift.size Numeric, the shift size for the moving window across the genes.
#' @param length.type Descriptive string of the length type being analyzed, typically "Gene".
#' @param methyl.type Descriptive string of the methylation type, typically something like "mCA/CA".
#' @param comp.between1 Descriptive string of the first condition or state being compared.
#' @param comp.between2 Descriptive string of the second condition or state being compared.
#'
#' @return A list containing ggplot objects for the overlay plots.
#' @export
#'
#' @examples
#' # Generating a toy dataset
#' set.seed(123) # For reproducibility
#' a <- runif(1000, min=-2, max=2)
#' b <- runif(1000, min=-2, max=2)
#' d <- runif(1000, min=0, max=0.05)
#' c <- sample(2000:1000000, 1000, replace=TRUE)
#' genes <- replicate(1000, paste(sample(letters, 8, replace = TRUE), collapse = ""))
#' df <- data.frame(comp.mat = a, logFC.crude = b, gene.length = c, mCA.CA = d, gene.name = genes)
#' # Running the overlaymC function
#' overlay_plots <- overlaymC(mat = df, bin.size = 60, shift.size = 6, methyl.type = "mCA/CA")
#' # To view the plots, you would call each element from the list, e.g., print(overlay_plots$plot1)

overlaymC <- function(mat, bin.size = 200, shift.size = 40,
                       length.type = "Gene", methyl.type = "",
                       comp.between1 = "", comp.between2 = ""){
  p1 <- overlaymCFunction(dat = mat, bin.size, shift.size,
                            length.type, methyl.type,
                            comp.between1, comp.between2)
  return(p1)
}

#' Generates Overlay Plots for Methylation and Gene Expression Analysis
#'
#' This function creates a series of overlay plots to explore the relationship 
#' between gene length, methylation values (mCA/CA), and gene expression levels.
#' It bins the genes according to length and calculates mean expression levels and 
#' methylation for each bin to identify trends across gene length categories.
#'
#' @param dat A dataframe where the first column is comp.mat representing 
#' comparison matrix values, the second column is log FC for log fold change, 
#' the third column is gene length, the fourth column contains the mCA.CA values, 
#' and the last column has the gene names.
#' @param bin.size Numeric, specifies the size of the bin for the analysis.
#' @param shift.size Numeric, indicates the shift size for the moving window 
#' in the binning process.
#' @param length.type A string indicating the type of length being analyzed 
#' (typically "Gene").
#' @param methyl.type A string indicating the type of methylation data being 
#' analyzed.
#' @param comp.between1 A string describing the first condition or group in 
#' the comparison.
#' @param comp.between2 A string describing the second condition or group in 
#' the comparison.
#'
#' @return A list of ggplot objects, each representing a different type of 
#' overlay plot for visualizing the binned methylation and gene expression data.
#' Each plot serves a unique purpose in the analysis:
#' - `overlapPlot`: Visualizes the relationship between mean log fold change 
#'   and methylation values across all bins.
#' - `mCAvGl_plot`: Shows how gene length varies with methylation values.
#' - `plotBar`: A bar plot indicating the distribution of gene lengths within each bin.
#' - `mCAvGl_inv_plot`: Similar to `mCAvGl_plot` but with inverted axes.
#' - `diffPlot`: Boxplot of the differences in methylation between long and short genes.
#' - `diff100kPlot`: Same as `diffPlot` but only for bins where the mean gene length is >= 100kb.
#' - `mat1`: A dataframe containing summary statistics for each bin.
#' - `mat2`: A list of dataframes, each corresponding to a bin and containing detailed gene data.
#' - `dat`: The original dataset after preprocessing.
#' @noRd
#' 
#' @examples
#' # Assuming that 'df' is a dataframe structured as required by the function:
#' result <- overlaymCFunction(dat = df, bin.size = 60, shift.size = 6, 
#'                             length.type = "Gene", methyl.type = "mCA/CA",
#'                             comp.between1 = "Condition1", comp.between2 = "Condition2")
#' # To view the first plot
#' print(result$overlapPlot)
#' # To access the summary statistics for each bin
#' head(result$mat1)


overlaymCFunction <- function(dat, bin.size, shift.size,
                                length.type, methyl.type,
                                comp.between1, comp.between2){
  dat <- dat[order(dat$mCA.CA),]
  #dat <- dat[order(dat$gene.length),]
  dat$gene.length <- dat$gene.length/1000
  dat$lg <- ifelse(dat$gene.length >= 100, "long gene", "short gene")
  dat$sublg <- ifelse(dat$gene.length >= 100, "long",
                      ifelse(dat$gene.length >= 50 & dat$gene.length < 100, "short",
                             "extra-short"))
  print(head(dat))
  ## binning the genes
  mean.points <- data.frame()
  mat.points <- list()
  num.bins <- round((dim(dat)[1]-bin.size)/shift.size)+1
  for(i in 0:num.bins){
    start <- i*shift.size+1
    end <- start + bin.size-1
    if ((start > dim(dat)[1])) break;
    if(end > dim(dat)[1]){
      end <- dim(dat)[1]
    }
    mat1 <- dat[start:end, 1]
    mat2 <- dat[start:end, 2]
    mat.points[[i+1]] <- dat[start:end, c(5,2:4,6)]
    mat.length <- mean(dat[start:end, 3])
    mC <- mean(dat[start:end, 4])
    mat.mean1 <- mean(mat1)
    mat.mean2 <- mean(mat2)
    mat.sd.1 <- sd(mat1)
    mat.sd.2 <- sd(mat2)
    bin.width <- end-start+1
    num.lg <- sum(dat[start:end, 3] >= 100)
    num.sm1 <- sum(dat[start:end, 3] >= 50 &
                     dat[start:end, 3] < 100)
    num.sm2 <- sum(dat[start:end, 3] < 50)

    ## including it to check for simpsons paradox
    dat0 <- dat[start:end, ]
    mat.mean.lg <- mean(as.numeric(dat0[which(dat0$sublg %in% "long"), 4]))
    mat.mean.sm1 <- mean(as.numeric(dat0[which(dat0$sublg %in% "short"), 4]))
    mat.mean.sm2 <- mean(as.numeric(dat0[which(dat0$sublg %in% "extra-short"), 4]))

    ## summary df
    mat.mean <- data.frame(mat.mean1, mat.mean2, mat.sd.1, mat.sd.2,
                           bin.width, mat.length, mC, num.lg, num.sm1, num.sm2,
                           mat.mean.lg, mat.mean.sm1, mat.mean.sm2)
    mean.points <- rbind(mean.points, mat.mean)
    if (end == dim(dat)[1]) break;
  }
  col1 <- comp.between1
  col2 <- comp.between2

  ## check the directionality for last 5 bins
  idx1 <- nrow(mean.points)-5
  idx2 <- nrow(mean.points)
  sign1 <- sign(mean.points$mat.mean1[idx1:idx2])
  sign2 <- sign(mean.points$mat.mean2[idx1:idx2])
  sign.comp <- sum(apply(data.frame(sign1, sign2), 1, function(r)(r[1] == r[2])))
  if(sign.comp <= 4){
    text1 <- "Warning: Directionality problem\nTherefore, inter-change "
    text2 <- "numerator and denominator for the calculation of logFC of control samples \n"
    message("\n",text1,text2)
    mean.points$mat.mean1 = -mean.points$mat.mean1
  }
  mean.points$pval <- apply(mean.points, 1, function(r){
    studentTest2(m1 = r[1], m2 = r[2],
            s1 = r[3], s2 = r[4], n1 = r[5])})
  mean.points$pval.log10 <- -log10(mean.points$pval)
  mean.points$fdr <- p.adjust(p = mean.points$pval, method = "fdr")
  mean.points$gene.type <- ifelse(test = mean.points$fdr < 0.05, 1, 0)
  mean.points$gl <- ifelse(##mean.points$mat.length >= 100 &
    mean.points$fdr < 0.05, "*",
    ifelse(mean.points$mat.length >= 100, "", ""))
  mean.points$gene.type <- factor(mean.points$gene.type, levels = c(1,0))

  ## print statements
  cat("Total number of bins = ", dim(mean.points)[1],"\n")
  cat("Total number of Short Gene bins = ",sum(mean.points$mat.length<100),"\n")
  cat("Total number of Long Gene bins = ",sum(mean.points$mat.length>=100),"\n")
  cat("Total number of bins that are statistically significant = ",
      sum(mean.points$gene.type == 1),"\n")
  cat("Total number of Short Gene bins that are statistically significant = ",
      sum(mean.points$gene.type == 1 & mean.points$mat.length < 100),"\n")
  cat("Total number of Long Gene bins that are statistically significant = ",
      sum(mean.points$gene.type == 1 & mean.points$mat.length >= 100),"\n")

  ## mC vs Log2FC
  plot1 <- ggplot(data = mean.points, aes(x = mC)) + theme_bw() +
    geom_line(aes(y = mat.mean1, color = col1), linewidth = 1) +
    geom_line(aes(y = mat.mean2, color = col2), linewidth = 1) +
    ylab(paste("Mean Log2FC")) +
    xlab(paste("Gene body", methyl.type)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50") +
    geom_ribbon(aes(ymin=(mat.mean2-(mat.sd.2*0.5)),
                    ymax=(mat.mean2+(mat.sd.2*0.5)),
                    x = mC, fill = col2), alpha=.25) +
    geom_ribbon(aes(ymin=(mat.mean1-(mat.sd.1*0.5)),
                    ymax=(mat.mean1+(mat.sd.1*0.5)),
                    x = mC, fill = col1), alpha=.25) +
    geom_point(aes(y = mat.mean2), shape=mean.points$gl, size=6) +
    theme(axis.title = element_text(size = 24, face = "bold",
                                    color = "black"),
          axis.text = element_text(size = 24, face = "bold",
                                   color = "black"),
          legend.position="none")

  ## mean gene length vs  mean mC
  idx1 <- which.max(mean.points[,"mat.length"])
  idx2 <- which.max(mean.points[,"num.lg"])
  cat("Bin with max mean gene length:",idx1,"\n")
  cat("Bin with max num of long genes (>= 100 Kb):",idx2,"\n")
  plot2 <- ggplot(data = mean.points, aes(y = mat.length, x = mC)) +
    geom_point() + theme_bw() + xlab(paste("Gene body", methyl.type)) +
    ylab(paste("Mean",length.type,"Length in KB")) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50")+
    geom_vline(xintercept = mean.points[idx1, "mC"],
               linetype = "longdash", color = "red") +
    geom_vline(xintercept = mean.points[idx2, "mC"],
               linetype = "dotted", color = "red") +
    annotate(geom="text", x=mean(mean.points[,"mC"]), y=50, color="black",
             label=paste("-- Max mean gene length Bin#",idx1,"\n",
                         ".. Max num of long genes Bin#",idx2),
             size=4, hjust=0, fontface = "bold")+
    theme(axis.title.y = element_text(size = 24, face = "bold",
                                      color = "black"),
          axis.text.y = element_text(size = 24, face = "bold",
                                     color = "black"),
          axis.title.x = element_text(size = 24, face = "bold",
                                      color = "black"),
          axis.text.x = element_text(size = 24, face = "bold",
                                     color = "black"),
          legend.position="none")

  ## barplot
  mean.points$bin <- paste0("Bin_",rownames(mean.points))
  mean.points$bin <- factor(mean.points$bin,
                            levels = stringi::stri_sort(mean.points$bin,
                                                        numeric = TRUE))
  melt.points <- melt(mean.points[, c("bin","num.lg","num.sm1","num.sm2")])
  plot3 <- ggplot(data=melt.points) +
    geom_bar(aes(y=value, x=bin, fill=variable), stat="identity") + ## position="stack",
    xlab("Bin Index") + ylab("Number of Genes") + theme_bw() +
    geom_hline(yintercept=bin.size/2, linetype="dashed", color="black")+
    geom_vline(xintercept=mean.points[idx1, "bin"], color="black",
               linetype="longdash") +
    geom_vline(xintercept=mean.points[idx2, "bin"], color="black",
               linetype="dotted") +
    annotate(geom="text",label=paste("Max mean gene length Bin#",idx1),
             x=(idx1-10), y=bin.size/2, angle=90, color="black",
             size=4) +
    annotate(geom="text",label=paste("Max num of long genes Bin#",idx2),
             x=(idx2+2), y=bin.size/2, angle=90, color="black",
             size=4) +
    scale_fill_discrete(name="Gene\nLength",
                        breaks=c("num.lg", "num.sm1", "num.sm2"),
                        labels=c(">= 100Kb", ">= 50kb & < 100kb",
                                 "< 50Kb")) +
    theme(axis.title.y = element_text(size = 24, face = "bold",
                                      color = "black"),
          axis.text.y = element_text(size = 24, face = "bold",
                                     color = "black"),
          axis.title.x = element_text(size = 24, face = "bold",
                                      color = "black"),
          legend.text = element_text(size = 14, face = "bold",
                                     color = "black"),
          legend.title = element_text(size = 14, face = "bold",
                                      color = "black"),
          axis.ticks.x = element_blank(), axis.text.x = element_blank())

  ## mC vs gene length -- inverted as compared to plot2
  plot4 <- ggplot(data = mean.points, aes(y = mC, x = mat.length)) +
    #scale_y_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
    geom_point() + theme_bw() + ylab(paste("Gene body", methyl.type)) +
    xlab(paste("Mean",length.type,"Length in KB")) +
    theme(axis.title.y = element_text(size = 24, face = "bold",
                                      color = "black"),
          axis.text.y = element_text(size = 24, face = "bold",
                                     color = "black"),
          axis.title.x = element_text(size = 24, face = "bold",
                                      color = "black"),
          axis.text.x = element_text(size = 24, face = "bold",
                                     color = "black"),
          legend.position="none")

  ## comparing bins where mean mCA is present for both long & short genes
  mean.points1 <- mean.points[,c("bin", "mC", "mat.mean.lg", "mat.mean.sm1",
                                 "mat.mean.sm2", "mat.length")]
  idx <- which(!is.na(mean.points1$mat.mean.lg) &
                 !is.na(mean.points1$mat.mean.sm1) &
                 !is.na(mean.points1$mat.mean.sm2))
  mean.points1 <- mean.points1[idx,]
  mean.points1$mat.mean.sm <- apply(mean.points1[,c("mat.mean.sm1",
                                                    "mat.mean.sm2")], 1, mean)
  mean.points1$mCA.diff <- apply(mean.points1[,c("mat.mean.lg", "mat.mean.sm")],
                                 1, function(r){r[1]-r[2]})
  plot5 <- ggplot(data = mean.points1, aes(y = mCA.diff, x = "")) +
    geom_boxplot() + geom_hline(yintercept = 0, linewidth = 0.75, alpha = 0.5,
                                linetype = "dashed", color = "red") +
    ylab("mCA diff (long - short) genes") + theme_bw() +
    theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
          axis.text.y = element_text(size = 24, face = "bold", color = "black"),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position = "none")

  plot6 <- ggplot(data = mean.points1[mean.points1$mat.length >= 100,],
                  aes(y = mCA.diff, x = "")) + geom_boxplot() +
    geom_hline(yintercept = 0, linewidth = 0.75, alpha = 0.5,
               linetype = "dashed", color = "red") +
    ylab("mCA diff (long - short)  genes") + theme_bw() +
    theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
          axis.text.y = element_text(size = 24, face = "bold", color = "black"),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position = "none")

  #colnames(mean.points) <- c()
  return(list(overlapPlot = plot1, mCAvGl_plot = plot2, plotBar = plot3,
              mCAvGl_inv_plot = plot4, diffPlot = plot5, diff100kPlot = plot6,
              mat1 = mean.points, mat2 = mat.points, dat = dat))
}

