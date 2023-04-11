## overlay plot
#' Overlay Gabel's plot
#'
#' @param mat data
#' @param bin.size bin size
#' @param shift.size shift size
#' @param comp.between1 comp between
#' @param comp.between2 comp between
#'
#' @return gabels plot
#' @export
#'
#' @examples
overlay.gabels.plot <- function(mat, bin.size = 200, shift.size = 40,
                                comp.between1 = "", comp.between2 = ""){

  p1 <- overlay.moving.average.function(dat = mat, bin.size, shift.size,
                                        comp.between1, comp.between2)
  return(p1)
}
#' Overlay Gabel's moving average function
#'
#' @param dat matrix
#' @param bin.size bin size
#' @param shift.size shift size
#' @param comp.between1 comp between
#' @param comp.between2 comp between
#'
#' @return overlay ggplots objects
#' @noRd
#'
#' @examples
overlay.moving.average.function <- function(dat, bin.size, shift.size,
                                            comp.between1, comp.between2){
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
  sign.comp <- sum(apply(data.frame(sign1, sign2), 1, function(r)(r[1] == r[2])))
  if(sign.comp < 3){
    text1 <- "Warning: Directionality issue\nTherefore, inter-change "
    text2 <- "numerator and denominator for the calculation of logFC of control samples \n"
    message("\n",text1,text2)
    mean.points$mat.mean1 = -mean.points$mat.mean1
  }

  ## calculating the p-value using studentt.test2
  mean.points$pval <- apply(mean.points, 1, function(r){
    studentt.test2(m1 = r[1], m2 = r[2],
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
    geom_ribbon(aes(ymin=(mat.mean1-(mat.sd.1*0.50)),
                    ymax=(mat.mean1+(mat.sd.1*0.50)),
                    x = mat.length, fill = col1), alpha=.25) +
    geom_ribbon(aes(ymin=(mat.mean2-(mat.sd.2*0.50)),
                    ymax=(mat.mean2+(mat.sd.2*0.50)),
                    x = mat.length, fill = col2), alpha=.25) +
    theme(## legend.text = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 24, face = "bold", color = "black"),
      axis.text.y = element_text(size = 24, face = "bold", color = "black"),
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      legend.position="none")

  ## output and p-value plot
  gene.type <- ifelse(test = mean.points$fdr < 0.05, yes = "#FF0000", no = "#696969")
  mean.points$gene.type <- ifelse(test = mean.points$fdr < 0.05, 1, 0)
  cat("Total number of bins = ", dim(mean.points)[1],"\n")
  cat("Total number of Short Gene bins = ", sum(mean.points$mat.length < 100),"\n")
  cat("Total number of Long Gene bins = ", sum(mean.points$mat.length >= 100),"\n")
  cat("Total number of bins that are statistically significant = ",
      sum(mean.points$gene.type == 1),"\n")
  cat("Total number of Short Gene bins that are statistically significant = ",
      sum(mean.points$gene.type == 1 & mean.points$mat.length < 100),"\n")
  cat("Total number of Long Gene bins that are statistically significant = ",
      sum(mean.points$gene.type == 1 & mean.points$mat.length >= 100),"\n")
  if(sum(mean.points$fdr < 0.05) > 0){
    y.int <- min(mean.points[which(mean.points$fdr < 0.05 &
                                     !is.infinite(mean.points$fdr)), "pval.log10"])
    plot2 <- ggplot(data = mean.points, aes(x = mat.length, y = pval.log10)) +
      geom_line(linewidth = 0.4, colour = "gray70") +
      geom_point(size = 2, color = gene.type) +
      geom_hline(aes(yintercept = y.int),
                 colour="#FF0000", linetype="dashed", size = 1) +
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
                 linetype="dashed", size = 1) + theme_bw() +
      theme(axis.title = element_text(size = 24, face = "bold"),
            axis.text.x = element_text(size = 24, face = "bold", color = "black"),
            axis.text.y = element_text(size = 24, face = "bold", color = "black"),
            legend.position="none")
  }
  return(list(plot1 = plot1, plot2 = plot2, bins.stat = mean.points, bins.info = mean.info))
}

## p-values from 2 sample t-test; code adapted from http://bit.ly/2eqeYyO
#' p-values from 2 sample t-test
#'
#' @param m1 the sample means
#' @param m2 the sample means
#' @param s1 the sample standard deviations
#' @param s2 the sample standard deviations
#' @param n1 the same sizes
#' @param n2 the same sizes
#' @param m0 the null value for the difference in means to be tested for. Default is 0
#' @param equal.variance whether or not to assume equal variance. Default is FALSE.
#'
#' @return p-values
#' @noRd
#'
#' @examples
studentt.test2 <- function(m1,m2,s1,s2,n1,n2=n1,m0=0,equal.variance=FALSE){
  if(equal.variance==FALSE){
    se <- sqrt( (s1^2/n1) + (s2^2/n2))
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
