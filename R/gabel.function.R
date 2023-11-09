#' Plot Moving Average for Comparison Between Groups
#'
#' @description
#' `gabelsPlot` applies Gabel's method or the moving average method to a given dataset
#' for the purpose of comparing two groups, typically knockout (KO) versus wild type (WT).
#' It calculates the moving average of a specified value, such as log2 fold change, across
#' a range of gene lengths and generates a plot to visualize trends in the data.
#'
#' @param mat A dataframe containing the variables for analysis: typically log fold change 
#' and gene length.
#' @param length.type A string describing the type of length being analyzed, usually "Gene".
#' This parameter is used to label the x-axis of the plot.
#' @param comp.between A string that describes the comparison being made, which will be 
#' used in the plot title to specify the groups being compared (e.g., "KO vs WT").
#' @param y.axis A string to be used as the label for the y-axis of the plot, typically 
#' "Mean Log2 Fold Change".
#'
#' @return A plot object created by the base `plot` function, representing the moving 
#' average of the specified value across gene lengths.
#' @export
#'
#' @examples
#' # Generate toy data for log fold change and gene length
#' set.seed(123) # For reproducibility
#' logFC.crude <- runif(1000, min=-2, max=2)
#' gene.length <- sample(2000:1000000, 1000, replace=TRUE)
#' df <- data.frame(logFC.crude = logFC.crude, gene.length = gene.length)
#'
#' # Plot the moving average for the generated data
#' gabelsPlot(mat = df, length.type = "Gene", 
#'            comp.between = "KO vs WT", 
#'            y.axis = "Mean Log2 Fold Change")

gabelsPlot <- function(mat, length.type = "Gene", comp.between = "",
                        y.axis = "Mean Log2 Fold Change"){
  p1 <- movingAverageFunction(dat = mat, bin.size = 200, shift.size = 40,
                                length.type, comp.between, y.axis)
  plot(p1)
}

#' Calculate and Plot Moving Average
#'
#' @description
#' Computes the moving average for a specified value across binned ranges of another variable, 
#' typically gene length. This is often used to analyze trends in gene expression data, such as 
#' log2 fold change across gene lengths. The function orders the data by gene length, bins it, 
#' calculates the mean for each bin, and then generates a ggplot object visualizing the trend.
#'
#' @param dat A dataframe or matrix where the first column contains the values for which the 
#' moving average should be calculated (e.g., log2 fold change) and the second column contains 
#' the variable used for binning (e.g., gene length).
#' @param bin.size The number of data points in each bin over which the moving average should 
#' be calculated.
#' @param shift.size The number of data points by which the bin window should shift for each 
#' iteration of the moving average calculation.
#' @param length.type The description of the binning variable, used for labeling the x-axis 
#' of the plot. Defaults to "Gene".
#' @param comp.between A descriptive label for the comparison being made, used in the y-axis 
#' label of the plot.
#' @param y.axis A descriptive label for the value being averaged, used in the y-axis label 
#' of the plot. Defaults to "Mean Log2 Fold Change".
#'
#' @return A ggplot object that displays the moving average plot.
#' @noRd
#'
#' @examples
#' # Assuming 'dat' is a dataframe with log2 fold change and gene length
#' example_dat <- data.frame(
#'   logFC = rnorm(1000, mean = 0, sd = 1),
#'   geneLength = runif(1000, min = 2000, max = 1000000)
#' )
#' moving_average_plot <- movingAverageFunction(
#'   dat = example_dat,
#'   bin.size = 50,
#'   shift.size = 10,
#'   length.type = "Gene",
#'   comp.between = "KO vs WT",
#'   y.axis = "Mean Log2 Fold Change"
#' )
#' print(moving_average_plot)

movingAverageFunction <- function(dat, bin.size, shift.size, length.type,
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
    mat.mean <- data.frame(mat.mean, mat.length)
    mean.points <- rbind(mean.points, mat.mean)
    ## end exceeds total number of genes
    if (end == dim(dat)[1]) break;
  }

  ## colors used for moving average plots
  col1 <- "#000000"
  col2 <- "#1E90FF"
  ind <- mean.points$mat.length >=1 & mean.points$mat.length <=1000
  mean.points <- mean.points[ind, ]
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
