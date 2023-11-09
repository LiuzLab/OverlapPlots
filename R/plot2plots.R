#' Normalization Check Box Plot
#'
#' @description
#' Generates a box plot to visually inspect whether the data is normalized. Each box represents the distribution
#' of data values for a sample, allowing for comparison across samples. Notches represent the confidence interval
#' around the median which can be used to assess the variance and skewness of the data distribution.
#'
#' @param data A numeric matrix or data frame where rows correspond to observations and columns to variables.
#' @param samples A factor vector with the same length as the number of columns in `data`.
#'                Each level corresponds to a group in the data that will be represented with a different color in the box plot.
#'
#' @return Invisibly returns an object of class `"boxplot"` which is a list with components that describe the box plot.
#'         The function also produces a box plot as a side effect.
#' @noRd
#'
#' @examples
#' # Assuming 'data' is a matrix of normalized values and 'samples' is a factor vector indicating groups
#' boxPlot(data = matrix(rnorm(100), ncol=10), samples = factor(rep(1:2, each=5)))
#'
boxPlot <- function(data, samples) {
  par(mar=c(2+round(max(nchar(colnames(data)))/2), 4, 2, 1), font=2)
  bp <- boxplot(data, boxwex=0.6, notch=TRUE, outline=FALSE, las=2, col=samples)
  legend("topleft", levels(samples), fill=palette(), bty="n", cex=0.8, xpd=TRUE)
  invisible(bp)
}

#' Gene Expression Level Visualization in Microarray Data
#'
#' @description
#' Generates a bar plot of gene expression levels for different samples and genotypes in microarray data using ggplot2.
#' It is designed to visualize the distribution of expression levels for a specific gene across various conditions.
#'
#' @param data A numeric matrix or data frame with genes in rows and samples in columns.
#' @param gene.id A character vector where the first element is the row name in `data` corresponding to the gene of interest, and the second element is used for the plot title.
#' @param genotypes A factor vector indicating the genotype for each sample.
#'
#' @return An object of class `ggplot` representing the gene expression level plot.
#' @noRd
#'
#' @examples
#' # Assuming 'data' is a data frame of expression levels with row names as gene IDs and 'genotypes' is a factor vector of genotypes
#' GeneExpLevels(data = your_data, gene.id = c("GeneX", "Gene X"), genotypes = your_genotypes)
#'
GeneExpLevels <- function(data, gene.id, genotypes) {
  ## MeCP2 expression level
  ind <- which(rownames(data) == gene.id[1])
  if (length(ind) == 0) stop("Gene ID not found in the row names of the data.")
  
  matGene <- as.vector(as.matrix(data[ind,]))
  plotDat <- data.frame(sampleName = colnames(data),
                        normCounts = matGene,
                        genotype = genotypes)
  
  ggplotObject <- ggplot(plotDat, aes(x = sampleName, y = normCounts, fill = genotype)) +
                  geom_bar(stat="identity") +
                  ylab("Log2 Normalized Counts") +
                  ggtitle(paste("Barplot for Gene ", gene.id[2], sep="")) +
                  theme_bw() +
                  facet_grid(. ~ genotype, space = "free", scales = "free") +
                  theme(plot.title = element_text(size = 24, face = "bold"),
                        axis.title.y = element_text(size = 20, color = "black", face = "bold"),
                        axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                        axis.text.x = element_blank(),
                        legend.position = "none",
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        strip.text = element_text(size = 18, color = "red", face = "bold"))
  
  return(ggplotObject)
}


#' Multidimensional Scaling (MDS) Plot
#'
#' @description
#' Creates a two-dimensional MDS plot using ggplot2. MDS is used to visualize the level of similarity of individual cases of a dataset.
#' The function scales the data and plots it, coloring and shaping points by genotypes and conditions.
#'
#' @param data A numeric matrix or data frame with rows as samples and columns as features.
#' @param genotypes A factor vector indicating the genotype for each sample, which will be used to color the points on the plot.
#' @param conditions An optional factor vector indicating the condition for each sample, which will be used to shape the points on the plot.
#'
#' @return A `ggplot` object representing the MDS plot.
#' @noRd
#'
#' @examples
#' # Assuming 'data' is a data frame with rows as genes and columns as samples
#' # 'genotypes' and 'conditions' are factor vectors indicating the genotype and condition for each sample, respectively
#' MDSplot(data = your_data, genotypes = your_genotypes, conditions = your_conditions)
#'
MDSplot <- function(data, genotypes, conditions){
  # Perform MDS
  mdsDist <- cmdscale(d = dist(t(data)), eig = TRUE, k = 2)
  mdsDist <- data.frame(genotypes, x = mdsDist$points[,1]/1e4,
                       y = mdsDist$points[,2]/1e4)

  # Create the base plot object
  base_plot <- ggplot(mdsDist, aes(x = x, y = y, color = genotypes))

  # Add points to the plot; use different shapes if 'conditions' is provided
  if(missing(conditions)){
    plot <- base_plot +
            geom_point(size = 1) +
            ylab("MDS Coordinate 2 (x 1e4)") +
            xlab("MDS Coordinate 1 (x 1e4)")
  }
  else{
    plot <- base_plot +
            aes(shape = conditions) +
            geom_point(size = 8) +
            ylab("MDS Coordinate 2 (x 1e4)") +
            xlab("MDS Coordinate 1 (x 1e4)")
  }

  # Customize the theme of the plot
  plot + theme_bw() +
         theme(legend.text = element_text(size = 18, face = "bold"),
               legend.title = element_text(size = 18, color = "black", face = "bold"),
               axis.title = element_text(size = 18, face = "bold"),
               axis.text.x = element_text(size = 18, face = "bold", color = "black"),
               axis.text.y = element_text(size = 18, face = "bold", color = "black"),
               plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
}


## PCA plot
#' Principal Component Analysis (PCA) Plot
#'
#' @description
#' Creates a PCA plot to visualize the first two principal components derived from the provided data.
#' PCA is used to reduce dimensionality and visualize the data in a two-dimensional space.
#'
#' @param data A numeric matrix with columns as samples and rows as features.
#' @param genotypes A factor that indicates the genotype of each sample, used to color the points in the plot.
#' @param conditions An optional factor indicating the condition of each sample, used to shape the points in the plot.
#'
#' @return A `ggplot` object representing the PCA plot with the first two principal components.
#' @noRd
#'
#' @examples
#' # Assuming 'data' is a numeric matrix with columns as samples
#' PCAplot(data = matrix(rnorm(100), ncol = 10), genotypes = factor(rep(1:2, each = 5)), conditions = factor(rep(1:2, each = 5)))
#'
PCAplot <- function(data, genotypes, conditions) {
  # Calculate the principal components
  pcs <- prcomp(t(data), center = TRUE)
  # Calculate the percentage of variance explained by the PCs
  percentVar <- round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2) * 100), 2)
  # Convert the PCA results into a data frame for plotting
  pcaData <- as.data.frame(pcs$x)
  
  # Define a helper function for making axis labels
  makeLab <- function(pv, pc) {
    paste0("PC", pc, " (", pv, "%)")
  }
  
  # Create the PCA plot
  if(missing(conditions)) {
    ggplot(pcaData, aes(PC1, PC2, color = genotypes)) +
      geom_point(size = 8) +
      xlab(makeLab(percentVar[1], 1)) + 
      ylab(makeLab(percentVar[2], 2))
  } else {
    ggplot(pcaData, aes(PC1, PC2, color = genotypes, shape = conditions)) +
      geom_point(size = 8) +
      xlab(makeLab(percentVar[1], 1)) + 
      ylab(makeLab(percentVar[2], 2))
  } + 
    theme_bw() +
    theme(legend.text = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 16, color = "black", face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold", color = "black"),
          axis.text.y = element_text(size = 16, face = "bold", color = "black"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
}


#' Create Labels for Principal Components with Variance
#'
#' @description
#' Generates a label string for principal components indicating the percentage of variance each component explains.
#'
#' @param x The percentage of variance explained by the principal component.
#' @param pc The principal component number.
#'
#' @return A character string labeling the principal component with its associated variance.
#' @noRd
#'
#' @examples
#' makeLab(25.3, 1)
#' # [1] "PC1: 25.3% variance"
#'
makeLab <- function(x, pc) {
  # Validate inputs
  stopifnot(is.numeric(x), x > 0, x <= 100, is.numeric(pc), pc > 0)
  
  # Create the label
  paste0("PC", pc, ": ", x, "% variance")
}

#' Scatter Plot Visualization for Differentially Expressed Genes (DEGs)
#'
#' @description
#' Generates a scatter plot of genes categorized by their expression level significance and length. 
#' Genes are classified as long or short based on gene length, and as up- or down-regulated based on log2 fold change.
#'
#' @param dat A data frame containing gene expression data with columns for gene name, log2 fold change, adjusted p-value, and gene length.
#' @param log2FC Numeric threshold for log2 fold change to consider a gene significantly up- or down-regulated.
#' @param comp.between A character string indicating the comparison between conditions or time points.
#' @param pval The p-value threshold to determine the statistical significance of gene expression changes.
#'
#' @return A ggplot object representing the scatter plot with annotations for the number of significantly up- and down-regulated long and short genes.
#' @noRd
#'
#' @examples
#' # Assuming 'dat' is a data frame with the necessary structure:
#' plotScatter(dat = your_data, log2FC = 1, comp.between = "Condition1 vs Condition2")
#'
plotScatter <- function(dat, log2FC, comp.between, pval = 0.05){
  # Setting up the data by renaming columns and categorizing genes
  colnames(dat) <- c("gene.name", "logFC", "adj.P.Val", "gene.length")
  gene.type <- ifelse(dat$adj.P.Val < pval & abs(dat$logFC) > log2FC,
                      ifelse(dat$gene.length > 100e3, "Long Genes", "Short Genes"),
                      "Not Stat. Signif.")
  ind <- which(gene.type != "Not Stat. Signif.")
  dat <- dat[ind,]
  gene.type <- gene.type[ind]

  # Generating a contingency table for the annotations
  up.LongGenes <- sum(dat$logFC > log2FC & dat$gene.length > 100e3)
  down.LongGenes <- sum(dat$logFC < -log2FC & dat$gene.length > 100e3)
  up.ShortGenes <- sum(dat$logFC > log2FC & dat$gene.length <= 100e3)
  down.ShortGenes <- sum(dat$logFC < -log2FC & dat$gene.length <= 100e3)
  cont.tab <- matrix(data = c(up.LongGenes, down.LongGenes, up.ShortGenes, down.ShortGenes),
                     nrow = 2, byrow = TRUE, dimnames = list(c("Up", "Down"), c("Long Genes", "Short Genes")))
  
  # Creating the scatter plot
  ggplot(data = dat, aes(y = logFC, x = gene.length/1000, colour = gene.type)) +
    geom_point(size = 1) +
    xlab("Gene Length (KB)") +
    ylab(paste("Log2 Fold Change", comp.between)) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,10,100,1000)) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    theme_bw() +
    annotate("text", x = 500, y = 1.5, label = cont.tab[1, 1], size = 5, fontface = "bold") +
    annotate("text", x = 500, y = -1.5, label = cont.tab[2, 1], size = 5, fontface = "bold") +
    annotate("text", x = 1, y = 1.5, label = cont.tab[1, 2], size = 5, fontface="bold") +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
            legend.position="none",
            axis.text.x = element_text(size = 18, face = "bold",
                                       color = "black"),
            axis.text.y = element_text(size = 18, face = "bold",
                                       color = "black"),
            legend.text = element_text(size = 18, face = "bold"),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
}


#' Scatter Plot with Linear Model Line and R-Squared Annotation
#'
#' @description
#' Produces a scatter plot with a linear regression line depicting the relationship
#' between gene length and log2 fold change. It calculates the R-squared value for the linear
#' model and annotates it on the plot.
#'
#' @param dat A data frame containing at least two numeric columns: `gene.length` and `logFC.crude`.
#'
#' @return A ggplot object showing the scatter plot with a linear model line and R-squared annotation.
#' @noRd
#'
#' @examples
#' # Assuming 'dat' is a data frame with columns 'gene.length' and 'logFC.crude':
#' scatterlm(dat = your_data)
#'
scatterlm <- function(dat){
  # Validate that necessary columns exist
  requiredCols <- c("gene.length", "logFC.crude")
  if(!all(requiredCols %in% names(dat))) {
    stop("Data frame must contain columns 'gene.length' and 'logFC.crude'")
  }
  
  # Fit linear model
  lm_fit <- lm(gene.length ~ logFC.crude, dat)
  
  # Calculate R-squared value and create label
  r.sq <- paste("R^2 = ", format(summary(lm_fit)$r.squared, digits = 2))
  
  # Print summary to console
  print(summary(lm_fit))
  
  # Create scatter plot with regression line
  scatter_plot <- ggplot(dat, aes(x = gene.length/1000, y = logFC.crude)) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue") +
    xlab("Gene Length (KB)") +
    ylab("Mean Log2 Fold Change") +
    theme_bw() +
    scale_x_continuous(trans = 'log10', breaks = c(0, 1, 10, 100, 1000)) +
    annotate("text", x = 500, y = 1.5, label = r.sq, size = 5, fontface = "bold") +
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
          axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
          axis.text.y = element_text(size = 18, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  # Return the plot
  return(scatter_plot)
}
