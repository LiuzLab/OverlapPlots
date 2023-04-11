## gabel's method or moving average method
#' Gabel Function
#'
#' @description
#' Gabel's method or moving average method
#'
#' @param mat data
#' @param length.type Gene
#' @param comp.between comp between
#' @param y.axis Mean Log2 Fold Change
#'
#' @return Box plot
#' @export
#'
#' @examples
gabels.plot <- function(mat, length.type = "Gene", comp.between = "",
                        y.axis = "Mean Log2 Fold Change"){
  p1 <- moving.average.function(dat = mat, bin.size = 200, shift.size = 40,
                                length.type, comp.between, y.axis)
  plot(p1)
}

#' Moving Average Function
#'
#' @param dat data
#' @param bin.size bin size
#' @param shift.size shift size
#' @param length.type Gene
#' @param comp.between comp between
#' @param y.axis Mean Lof2 Fold Change
#'
#' @return ggplot object
#' @noRd
#'
#' @examples
moving.average.function <- function(dat, bin.size, shift.size, length.type,
                                    comp.between, y.axis){
  dat <- dat[order(dat[,2]),]
  dat[,2] <- dat[,2]/1000

  # data frame for storing the values and
  # calculating the number of bins
  mean.points <- data.frame()
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
    mat <- dat[start:end, 1]
    mat.mean <- mean(mat)
    mat.length <- mean(dat[start:end, 2])
    #cat("Bin",i+1,": ",start,"-",end,"\t",mat.length,"\t",mat.mean,"\n",sep="")
    mat.mean <- data.frame(mat.mean, mat.length)
    mean.points <- rbind(mean.points, mat.mean)
    ## end exceeds total number of genes
    if (end == dim(dat)[1]) break;
  }

  ## colors used for moving average plots
  col1 <- "#000000"
  col2 <- "#1E90FF"
  ind = mean.points$mat.length >=1 & mean.points$mat.length <=1000
  mean.points = mean.points[ind, ]
  plot1 <- ggplot(data = mean.points, aes(x = mat.length, y = mat.mean)) +
    geom_point(size = 1.5, colour = col2) + geom_line(size=1, color = col1) +
    scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
    xlab(paste("Mean",length.type,"Length in KB")) +
    ylab(paste(y.axis, comp.between)) + theme_bw() +
    theme(legend.position = "none", plot.title = element_blank(),
          axis.title = element_text(size = 24, face = "bold"),
          axis.text.x = element_text(size = 22, face = "bold", color = "black"),
          axis.text.y = element_text(size = 22, face = "bold", color = "black"),
          plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
  return(plot1)
}
