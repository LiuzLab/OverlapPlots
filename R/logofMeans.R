#' Calculate Log2 Fold Change Between Two Sets of Samples
#'
#' @description
#' The function calculates the average expression values for two sets of samples (A and B) and 
#' then computes the fold change and log2 fold change between them. This is typically used 
#' to analyze the differential expression of genes between two conditions in RNA-seq data.
#'
#' @param dat A dataframe or matrix where rows represent genes and columns represent samples. 
#' The differential expression results from `DESeqCalculation` can be used as input.
#' @param A.samples A vector of column indices or names for the first sample type (e.g., control group).
#' @param B.samples A vector of column indices or names for the second sample type (e.g., treatment group).
#'
#' @return A dataframe with additional columns for the mean expression of samples A and B, 
#' crude fold change, and log2 fold change. Rows with `NA` values in the log2 fold change 
#' are removed from the output.
#' @export
#'
#' @examples
#' # Simulate RNA-seq count data for 100 genes and 6 samples
#' set.seed(123) # for reproducibility
#' dat <- as.data.frame(matrix(rnbinom(n=600, mu=100, size=0.5), ncol=6))
#' colnames(dat) <- c(paste("Sample_A", 1:3, sep=""), paste("Sample_B", 1:3, sep=""))
#' rownames(dat) <- paste("Gene", 1:100, sep="")
#'
#' # Calculate log2 fold change between the sample groups
#' logFC_results <- logofMeansBetweenAB(dat = dat, A.samples = 1:3, B.samples = 4:6)
#' head(logFC_results)


logofMeansBetweenAB <- function(dat, A.samples, B.samples){
  dat$Mean.A <- apply(dat[,A.samples], 1, function(r) {(mean(r))})
  dat$Mean.B <- apply(dat[,B.samples], 1, function(r) {(mean(r))})
  dat$FC.crude <- apply(dat[,c("Mean.A", "Mean.B")], 1,
                        function(r) {(r[2]/r[1])})
  dat$logFC.crude <- apply(dat[,c("Mean.A", "Mean.B")], 1,
                        function(r) {log2((r[2]+1)/(r[1]+1))})
  dat <- dat[!is.na(dat$logFC.crude),]
  return(dat)
}

#' Calculate Log2 Fold Change Among Three Sets of Samples
#'
#' @description
#' The function computes the average expression values for three sets of samples (A, B, and C) 
#' and then calculates the fold change and log2 fold change using the differences between 
#' these averages. This function is useful for RNA-seq data analysis where comparisons are 
#' made across three different conditions or time points.
#'
#' @param dat A dataframe or matrix with rows representing genes and columns representing 
#' samples. It typically contains count data from RNA-seq experiments.
#' @param A.samples A vector of column indices or names for the first sample group.
#' @param B.samples A vector of column indices or names for the second sample group.
#' @param C.samples A vector of column indices or names for the third sample group.
#'
#' @return A modified dataframe that includes additional columns for the mean expression 
#' of each sample group, the crude fold change, and the log2 fold change based on the 
#' differences between the sample means. Rows containing `NA` values after computation 
#' are excluded from the output.
#' @export
#'
#' @examples
#' # Simulate RNA-seq count data for 100 genes across 9 samples
#' set.seed(123) # for reproducibility
#' dat <- matrix(rnbinom(n=900, mu=100, size=0.5), ncol=9)
#' colnames(dat) <- c(paste("Sample_A", 1:3, sep=""), 
#'                    paste("Sample_B", 1:3, sep=""), 
#'                    paste("Sample_C", 1:3, sep=""))
#' rownames(dat) <- paste("Gene", 1:100, sep="")
#' dat <- as.data.frame(dat)
#'
#' # Calculate log2 fold change among the three sample groups
#' logFC_results <- logofMeansBetweenABC(dat = dat, A.samples = 1:3, 
#'                                       B.samples = 4:6, C.samples = 7:9)
#' head(logFC_results)

logofMeansBetweenABC <- function(dat, A.samples, B.samples, C.samples){
  dat$Mean.A <- apply(dat[,A.samples], 1, function(r) {(mean(r))})
  dat$Mean.B <- apply(dat[,B.samples], 1, function(r) {(mean(r))})
  dat$Mean.C <- apply(dat[,C.samples], 1, function(r) {(mean(r))})
  dat$FC.crude <- apply(dat[,c("Mean.A", "Mean.B", "Mean.C")], 1,
                        function(r) {((r[2]-r[1])/(r[3]-r[1]))})
  dat$logFC.crude <- apply(dat[,c("Mean.A", "Mean.B", "Mean.C")], 1,
                        function(r) {log2((r[2]-r[1]+1)/(r[3]-r[1]+1))})
  dat <- dat[!is.na(dat$logFC.crude),]
  return(dat)
}

#' Calculate Log2 Fold Changes Within Genotypes
#'
#' @description
#' Computes log2 fold changes (log2FC) between pairs of genotypes within the provided data. 
#' It first converts log2 values back to intensity values, calculates the mean intensities 
#' for each genotype, and then derives the log2FC between these means. This function is 
#' particularly useful for analyzing differential expression within genotypic comparisons.
#'
#' @param dat A dataframe where rows are genes and columns are samples, which must include 
#' 'gene.name' and 'gene.length' columns along with log2-transformed expression values for each genotype.
#'
#' @return A dataframe containing log2 fold change values for each pair of genotypes 
#' along with gene length information.
#'
#' @noRd
#'
#' @examples
#' # Simulated example with gene names and gene lengths
#' dat <- data.frame(gene.name = paste("Gene", 1:100, sep=""),
#'                   gene.length = sample(1000:2000, 100, replace = TRUE),
#'                   replicate(matrix(rnorm(200, mean = 5, sd = 0.3), ncol = 2), 2))
#' colnames(dat)[3:6] <- c("Genotype1_Rep1", "Genotype1_Rep2", 
#'                         "Genotype2_Rep1", "Genotype2_Rep2")
#'
#' # Calculate log2 fold change within genotypes
#' log2FC_results <- log2FCwithingenotypes(dat)
#' head(log2FC_results)

log2FCwithingenotypes <- function(dat){

  ## converting log2 values to intensity and then calculating log2FC
  log2FC.dat <- data.frame(row.names = rownames(dat))
  rownames(dat) <- dat$gene.name
  gene.length <- dat[,c("gene.length","gene.name")]
  dat <- dat[,!names(dat) %in% c("gene.name","gene.length")]
  for(j in 2:ncol(dat)){
    i = j-1
    if(j > i){
      sample1 <- colnames(dat[,c(i,j)])
      sample2 <- colnames(dat[,-c(i,j)])
      dat.mean <- logofMeansBetweenAB(dat = dat, A.samples = sample2,
                                         B.samples = sample1)
      log2FC.dat <- cbind(log2FC.dat, dat.mean$logFC.crude)
    }
  }
  colnames(log2FC.dat) <- c("comp.mat1", "comp.mat2", "comp.mat3")
  log2FC.dat <- data.frame(log2FC.dat, gene.length)
  return(log2FC.dat)
}
