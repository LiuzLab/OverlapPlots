#' Overlay Gabel's Plot for Gene Expression Analysis
#'
#' @description
#' Creates an overlay plot that displays a moving average of log fold change (log FC) 
#' across binned gene lengths. The plot includes a confidence interval around the moving average.
#' This type of visualization helps in identifying trends in gene expression related to gene length.
#'
#' @param mat A dataframe where the first column is `comp.mat` representing comparison 
#' matrix values, the second column is `log FC` for log fold change, and the third 
#' column is `gene.length`.
#' @param bin.size Numeric, specifies the size of the bin for the analysis.
#' @param shift.size Numeric, indicates the shift size for the moving window 
#' in the binning process.
#' @param comp.between1 String, description of the first condition or state being compared.
#' @param comp.between2 String, description of the second condition or state being compared.
#' @param confidenceinterval Numeric, the width of the confidence interval to be displayed 
#' on the plot.
#'
#' @return A ggplot object that represents the Gabel's plot with a moving average 
#' and confidence interval for gene expression analysis.
#' @export
#'
#' @examples
#' # Generating a toy dataset
#' set.seed(123) # For reproducibility
#' a <- runif(1000, min=-2, max=2)
#' b <- runif(1000, min=-2, max=2)
#' c <- sample(2000:1000000, 1000, replace=TRUE)
#' df <- data.frame(comp.mat = a, logFC.crude = b, gene.length = c)
#' # Running the overlayGabelsPlot function
#' gabels_plot <- overlayGabelsPlot(mat = df, comp.between1 = "(WT/WT)",
#'                                  comp.between2 = "(KO/WT)", bin.size = 200,
#'                                  shift.size = 40, confidenceinterval=0.50)
#' # To display the plot
#' print(gabels_plot)
overlayGabelsPlot <- function(mat, bin.size = 200, shift.size = 40,
                              comp.between1 = "", comp.between2 = "",
                              confidenceinterval = 0.50){

  p1 <- overlayMovingAverageFunction(dat = mat, bin.size, shift.size,
                              comp.between1, comp.between2, confidenceinterval)
  return(p1)
}

#' Overlay Moving Average Plot for Gene Expression Analysis
#'
#' @description
#' Creates an overlay plot displaying moving averages of gene expression levels (log FC) 
#' across binned gene lengths. The function orders genes by length, bins them, and calculates
#' the mean log fold change and standard deviation within each bin. The overlay plot includes
#' a line for the moving average and a shaded confidence interval around the average.
#'
#' @param dat A matrix where the first column is `comp.mat`, representing comparison matrix 
#' values, the second column is `log FC` for log fold change, and the third column is 
#' `gene.length`.
#' @param bin.size Numeric, the size of the bin for aggregating gene data.
#' @param shift.size Numeric, the shift size for the moving window across the gene length.
#' @param comp.between1 String, description of the first condition or group being compared.
#' @param comp.between2 String, description of the second condition or group being compared.
#' @param confidenceinterval Numeric, the confidence interval width for the moving average 
#' plot. A value of 0.50 corresponds to a 50% confidence interval.
#'
#' @return A list of ggplot objects that represent different aspects of the gene expression 
#' analysis. The `plot1` object shows the moving average lines with confidence intervals, 
#' and `plot2` displays the -log10(p-value) across the gene length bins. The list also 
#' contains a dataframe `bins.stat` with statistics for each bin and `bins.info` with 
#' information such as the start and end of each bin and the count of genes with positive 
#' and negative log fold changes.
#' @noRd
#'
#' @examples
#' # Assuming 'df' is a dataframe with the necessary structure:
#' results <- overlayMovingAverageFunction(dat = df, bin.size = 60, shift.size = 6,
#'                                         comp.between1 = "Condition1", 
#'                                         comp.between2 = "Condition2", 
#'                                         confidenceinterval = 0.50)
#' # To plot the moving average with confidence interval:
#' print(results$plot1)
#' # To plot the -log10(p-value) across gene length bins:
#' print(results$plot2)

overlayMovingAverageFunction <- function(dat, bin.size, shift.size,
                                         comp.between1, comp.between2,
                                         confidenceinterval){
  dat <- dat[order(dat$gene.length),]
  dat$gene.length <- dat$gene.length/1000

  # data frame for storing the values and
  # calculating the number of bins
  mean.points <- data.frame()
  mean.info <- data.frame()
  num.bins <- round((dim(dat)[1]-bin.size)/shift.size)+1

  ## taking the mean of log2FC and genomic length
  for(i in 0:num.bins){
    start <- i*shift.size+1
    end <- start + bin.size-1
    ## if the start exceeds total number of genes
    if ((start > dim(dat)[1])) break;

    ## if the last bin exceeds the number of genes available
    if(end > dim(dat)[1]){
      end <- dim(dat)[1]
    }
    mat1 <- dat[start:end, 1]
    mat2 <- dat[start:end, 2]
    mat.mean1 <- mean(mat1)
    mat.mean2 <- mean(mat2)
    mat.sd.1 <- sd(mat1)
    mat.sd.2 <- sd(mat2)
    mat.length <- mean(dat[start:end, 3])
    bin.width <- end-start+1
    pos.logfc <- sum(dat[start:end, 1] >= 0)
    neg.logfc <- sum(dat[start:end, 1] < 0)
    overall.logfc <- ifelse(pos.logfc >= neg.logfc, "yes", "no")

    ## mat means
    mat.mean <- data.frame(mat.mean1, mat.mean2, mat.sd.1, mat.sd.2,
                           bin.width, mat.length)
    mat.info <- data.frame(start, end, pos.logfc, neg.logfc, overall.logfc,
                           mat.length)
    mean.points <- rbind(mean.points, mat.mean)
    mean.info <- rbind(mean.info, mat.info)

    ## end exceeds total number of genes
    if (end == dim(dat)[1]) break;
  }
  ## colors used for moving average plots
  col1 <- comp.between1
  col2 <- comp.between2

  ## check the directionality for last 5 bins
  idx1 <- nrow(mean.points)-4
  idx2 <- nrow(mean.points)
  sign1 <- sign(mean.points$mat.mean1[idx1:idx2])
  sign2 <- sign(mean.points$mat.mean2[idx1:idx2])
  sign.comp <- sum(apply(data.frame(sign1, sign2), 1,
                         function(r)(r[1] == r[2])))
  if(sign.comp < 3){
    text1 <- "Warning: Directionality issue\nTherefore, inter-change "
    text2 <- paste("numerator and denominator for the calculation of logFC",
    "of control samples \n")
    message("\n",text1,text2)
    mean.points$mat.mean1 = -mean.points$mat.mean1
  }

  ## calculating the p-value using studentTest2

  mean.points$pval <- apply(mean.points, 1, function(r){
    studentTest2(m1 = r[1], m2 = r[2],
            s1 = r[3], s2 = r[4], n1 = r[5])})
  mean.points$pval.log10 <- -log10(mean.points$pval)

  ## overlay line plot
  mean.points$fdr <- p.adjust(p = mean.points$pval, method = "fdr")
  # ind <- mean.points$mat.length >=1 & mean.points$mat.length <=1000
  # mean.points = mean.points[ind, ]
  plot1 <- ggplot(data = mean.points, aes(x = mat.length)) +
    geom_line(aes(y = mean.points$mat.mean1, color = col1),
              linewidth = 1) +
    geom_line(aes(y = mean.points$mat.mean2, color = col2),
              linewidth = 1) +
    ylab(paste("Mean Log2FC")) + theme_bw() +
    scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
    geom_ribbon(aes(ymin=(mat.mean1-(mat.sd.1*confidenceinterval)),
                    ymax=(mat.mean1+(mat.sd.1*confidenceinterval)),
                    x = mat.length, fill = col1), alpha=.25) +
    geom_ribbon(aes(ymin=(mat.mean2-(mat.sd.2*confidenceinterval)),
                    ymax=(mat.mean2+(mat.sd.2*confidenceinterval)),
                    x = mat.length, fill = col2), alpha=.25) +
    theme(## legend.text = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 24, face = "bold", color = "black"),
      axis.text.y = element_text(size = 24, face = "bold", color = "black"),
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      legend.position="none")

  ## output and p-value plot
  gene.type <- ifelse(test = mean.points$fdr < 0.05, yes = "#FF0000",
                      no = "#696969")
  mean.points$gene.type <- ifelse(test = mean.points$fdr < 0.05, 1, 0)
  # Commenting cat functions out to avoid BiocConductor NOTES 4/14/2023 - DP
  # cat("Total number of bins = ", dim(mean.points)[1],"\n")
  # cat("Total number of Short Gene bins = ", sum(mean.points$mat.length < 100),
  #     "\n")
  # cat("Total number of Long Gene bins = ", sum(mean.points$mat.length >= 100),
  #     "\n")
  # cat("Total number of bins that are statistically significant = ",
  #     sum(mean.points$gene.type == 1),"\n")
  # cat("Total number of Short Gene bins that are statistically significant = ",
  #     sum(mean.points$gene.type == 1 & mean.points$mat.length < 100),"\n")
  # cat("Total number of Long Gene bins that are statistically significant = ",
  #     sum(mean.points$gene.type == 1 & mean.points$mat.length >= 100),"\n")
  if(sum(mean.points$fdr < 0.05) > 0){
    y.int <- min(mean.points[which(mean.points$fdr < 0.05 &
                                     !is.infinite(mean.points$fdr)),
                             "pval.log10"])
    plot2 <- ggplot(data = mean.points, aes(x = mat.length, y = pval.log10)) +
      geom_line(linewidth = 0.4, colour = "gray70") +
      geom_point(size = 2, color = gene.type) +
      geom_hline(aes(yintercept = y.int),
                 colour="#FF0000", linetype="dashed", linewidth = 1) +
      scale_x_continuous(trans = log10_trans(),
                         breaks = c(0,1,10,100,1000)) +
      xlab(paste("Mean Gene Length in KB")) +
      ylab(paste("-Log10(pvalue)")) + theme_bw() +
      theme(legend.position="none",
            axis.title = element_text(size = 24, face = "bold"),
            axis.text.x = element_text(size = 24, face = "bold",
                                       color = "black"),
            axis.text.y = element_text(size = 24, face = "bold",
                                       color = "black"))
  }else{
    y.int <- ceiling(max(mean.points[,"pval.log10"]))
    plot2 <- ggplot(data = mean.points, aes(x = mat.length, y = pval.log10)) +
      geom_line(linewidth = 0.4, colour = "gray70") +
      geom_point(size = 2, color = gene.type) +
      scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
      xlab(paste("Mean Gene Length in KB")) + ylab(paste("-Log10(pvalue)")) +
      geom_hline(aes(yintercept = y.int), colour="#FF0000",
                 linetype="dashed", linewidth = 1) + theme_bw() +
      theme(axis.title = element_text(size = 24, face = "bold"),
            axis.text.x = element_text(size = 24, face = "bold",
                                       color = "black"),
            axis.text.y = element_text(size = 24, face = "bold",
                                       color = "black"),
            legend.position="none")
  }
  return(list(plot1 = plot1, plot2 = plot2, bins.stat = mean.points,
              bins.info = mean.info))
}

#' Calculate P-Values from a Two-Sample T-Test
#'
#' This function performs a two-sample t-test to compare the means of two groups.
#' It calculates the t-statistic and the corresponding p-value, allowing for unequal 
#' variances between the groups by default, a scenario often referred to as Welch's t-test.
#'
#' @param m1 The mean of the first sample.
#' @param m2 The mean of the second sample.
#' @param s1 The standard deviation of the first sample.
#' @param s2 The standard deviation of the second sample.
#' @param n1 The size of the first sample.
#' @param n2 The size of the second sample, defaulting to the size of the first sample (`n1`) if not specified.
#' @param m0 The hypothesized difference between the two means under the null hypothesis, default is 0.
#' @param equal.variance A logical flag indicating whether to assume equal variances for the two samples, 
#' default is `FALSE`.
#'
#' @return The p-value resulting from the two-sample t-test.
#' @export
#'
#' @examples
#' # Generate random data for two samples
#' a <- runif(30, min=-0.2, max=0.2)
#' b <- runif(30, min=-0.2, max=0.2)
#' a_std <- runif(30, min=1, max=1.2)
#' b_std <- runif(30, min=1.3, max=1.5)
#' n <- rep(30, 30)
#' # Calculate p-values using the studentTest2 function
#' p_values <- mapply(studentTest2, m1 = a, m2 = b, s1 = a_std, s2 = b_std, n1 = n)
#' print(p_values)

studentTest2 <- function(m1,m2,s1,s2,n1,n2=n1,m0=0,equal.variance=FALSE){
  if(equal.variance==FALSE){
    se <- sqrt(abs((s1^2/n1) + (s2^2/n2)))
    df <- ((s1^2/n1 + s2^2/n2)^2)/((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
  }
  else{
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2))
    df <- n1+n2-2
  }
  t <- (m1-m2-m0)/se
  pval <- 2*pt(-abs(t),df)
  return(pval)
}
