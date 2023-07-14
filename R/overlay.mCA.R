## overlay plot

#' mC Overlay Plot
#'
#' @param mat dataframe with first column corresponding to comp.mat, second
#' column to log FC, and third columns to gene length.
#' @param bin.size bin size
#' @param shift.size shift size
#' @param length.type length type
#' @param methyl.type methylation type
#' @param comp.between1 comp between
#' @param comp.between2 comp between
#'
#' @return Overlay plots mC
#' @export
#'
#' @examples
overlay.mC <- function(mat, bin.size = 200, shift.size = 40,
                       length.type = "Gene", methyl.type = "",
                       comp.between1 = "", comp.between2 = ""){
  p1 <- overlay.mC.function(dat = mat, bin.size, shift.size,
                            length.type, methyl.type,
                            comp.between1, comp.between2)
  return(p1)
}
overlay.mC.function <- function(dat, bin.size, shift.size,
                                length.type, methyl.type,
                                comp.between1, comp.between2){
  dat <- dat[order(dat$mCA.CA),]
  #dat <- dat[order(dat$gene.length),]
  dat$gene.length <- dat$gene.length/1000
  dat$lg <- ifelse(dat$gene.length >= 100, "long gene", "short gene")

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
    mat.mean <- data.frame(mat.mean1, mat.mean2, mat.sd.1, mat.sd.2,
                           bin.width, mat.length, mC, num.lg,
                           num.sm1, num.sm2)
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
  cat("Total number of Short Gene bins = ", sum(mean.points$mat.length < 100),"\n")
  cat("Total number of Long Gene bins = ", sum(mean.points$mat.length >= 100),"\n")
  cat("Total number of bins that are statistically significant = ",
      sum(mean.points$gene.type == 1),"\n")
  cat("Total number of Short Gene bins that are statistically significant = ",
      sum(mean.points$gene.type == 1 & mean.points$mat.length < 100),"\n")
  cat("Total number of Long Gene bins that are statistically significant = ",
      sum(mean.points$gene.type == 1 & mean.points$mat.length >= 100),"\n")

  ## plot1 -- mC vs Log2FC
  plot1 <- ggplot(data = mean.points, aes(x = mC)) + theme_bw() +
                  geom_line(aes(y = mat.mean1, color = col1), linewidth = 1) +
                  geom_line(aes(y = mat.mean2, color = col2), linewidth = 1) +
                  ylab(paste("Mean Log2FC")) + xlab("") +
                  geom_hline(yintercept = 0, linetype = "dashed",
                             color = "grey50") +
                  geom_ribbon(aes(ymin=(mat.mean2-(mat.sd.2*0.5)),
                                  ymax=(mat.mean2+(mat.sd.2*0.5)),
                                  x = mC, fill = col2), alpha=.25) +
                  geom_ribbon(aes(ymin=(mat.mean1-(mat.sd.1*0.5)),
                                  ymax=(mat.mean1+(mat.sd.1*0.5)),
                                  x = mC, fill = col1), alpha=.25) +
                  geom_point(aes(y = mat.mean2), shape = mean.points$gl, size = 6) +
                  theme(axis.title.y = element_text(size = 24, face = "bold",
                                                    color = "black"),
                        axis.text.y = element_text(size = 24, face = "bold",
                                                   color = "black"),
                        axis.title.x = element_blank(), axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(), legend.position="none")

  ## plot2 -- mC vs gene length
  plot2 <- ggplot(data = mean.points, aes(y = mat.length, x = mC, color = gene.type)) +
            #scale_y_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
            geom_point() + theme_bw() + xlab(paste("Gene body", methyl.type)) +
            ylab(paste("Mean",length.type,"Length in KB")) +
            geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
            theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.title.x = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.x = element_text(size = 24, face = "bold", color = "black"),
                  legend.position="none")
  plot3 <- plot1 + xlab(paste("Gene body", methyl.type)) +
            theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.title.x = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.x = element_text(size = 24, face = "bold", color = "black"),
                  legend.position="none")

  ## barplot
  mean.points$bin <- paste0("Bin_",rownames(mean.points))
  mean.points$bin <- factor(mean.points$bin,
                            levels = stringi::stri_sort(mean.points$bin,
                                                        numeric = TRUE))
  melt.points <- melt(mean.points[##mean.points$mat.length >= 100,
                                  c("bin","num.lg","num.sm1","num.sm2")])
  plot4 <- ggplot(data=melt.points, aes(fill=variable, y=value, x=bin)) +
            geom_bar(position="stack", stat="identity") + theme_bw() +
            xlab("Bin Index") + ylab("Number of Genes") +
            geom_hline(yintercept = bin.size/2, linetype = "dashed", color = "black") +
            # geom_line(aes(x = mean.points[which(mean.points$mC >= 0.01)[1],"bin"]),
            #           color = "black", size = 1) +
            geom_line(aes(x = mean.points[which(mean.points$mC >= 0.02)[1],"bin"]),
                      color = "black", size = 1) +
            geom_line(aes(x = mean.points[which(mean.points$mC >= 0.03)[1],"bin"]),
                      color = "black", size = 1) +
            # geom_line(aes(x = mean.points[which(mean.points$mC >= 0.04)[1],"bin"]),
            #           color = "black", size = 1) +
            # geom_line(aes(x = mean.points[which(mean.points$mC >= 0.05)[1],"bin"]),
            #           color = "black", size = 1) +
            scale_fill_discrete(name="Gene\nLength",
                                breaks=c("num.lg", "num.sm1", "num.sm2"),
                                labels=c(">= 100Kb", ">= 50kb", "< 50Kb")) +
            theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.title.x = element_text(size = 24, face = "bold", color = "black"),
                  legend.text = element_text(size = 18, face = "bold", color = "black"),
                  legend.title = element_text(size = 18, face = "bold", color = "black"),
                  axis.ticks.x = element_blank(), axis.text.x = element_blank())
  plot5 <- ggplot(data = mean.points, aes(y = mC, x = mat.length, color = gene.type)) +
            #scale_y_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
            geom_point() + theme_bw() + ylab(paste("Gene body", methyl.type)) +
            xlab(paste("Mean",length.type,"Length in KB")) +
            theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.y = element_text(size = 24, face = "bold", color = "black"),
                  axis.title.x = element_text(size = 24, face = "bold", color = "black"),
                  axis.text.x = element_text(size = 24, face = "bold", color = "black"),
                  legend.position="none")
  plot(plot5 / plot4)
  return(list(plot1 = plot1, plot2 = plot2, plot3 = plot3, plot4 = plot4, plot5 = plot5,
              mat1 = mean.points, mat2 = mat.points, dat = dat))
}
