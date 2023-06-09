% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logofMeans.R
\name{logofMeansBetweenAB}
\alias{logofMeansBetweenAB}
\title{Log of means between two samples}
\usage{
logofMeansBetweenAB(dat, A.samples, B.samples)
}
\arguments{
\item{dat}{a dataframe. Differential analysis results from "DESeqCalculation"}

\item{A.samples}{First sample type}

\item{B.samples}{Second sample type}
}
\value{
log2 fold change of means between two given samples
}
\description{
Computes the mean of each sample, fold change, and log fold change.
}
\examples{
dat <- matrix(rnorm(600, sd = 0.3), ncol = 6)
dat <- as.data.frame(dat)
logofMeansBetweenAB(dat = dat, A.samples = 1:3, B.samples = 4:6)

}
