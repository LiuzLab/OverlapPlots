% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlay.gabel.function.R
\name{studentTest2}
\alias{studentTest2}
\title{p-values from 2 sample t-test}
\usage{
studentTest2(m1, m2, s1, s2, n1, n2 = n1, m0 = 0, equal.variance = FALSE)
}
\arguments{
\item{m1}{the sample means}

\item{m2}{the sample means}

\item{s1}{the sample standard deviations}

\item{s2}{the sample standard deviations}

\item{n1}{the same sizes}

\item{n2}{the same sizes}

\item{m0}{the null value for the difference in means to be tested for.
Default is 0}

\item{equal.variance}{whether or not to assume equal variance.
Default is FALSE.}
}
\value{
p-values
}
\description{
p-values from 2 sample t-test
}
\examples{
a <- runif(30, min=-0.2, max=0.2)
b <- runif(30, min=-0.2, max=0.2)
astd <- runif(30, min=1, max=1.2)
bstd <- runif(30, min=1.3, max=1.5)
binwidth <- rep(c(200), each = 30)
length <- runif(30, min=100, max=900)
df <- data.frame(mat.mean1 = a, mat.mean2 = b, mat.sd.1 = astd,
mat.sd.2 = bstd, bin.width = binwidth, mat.length = length)
r <- as.matrix(df)
r$pval <- apply(r, 1, function(r){
studentTest2(m1 = r[1], m2 = r[2],
            s1 = r[3], s2 = r[4], n1 = r[5])})
}
