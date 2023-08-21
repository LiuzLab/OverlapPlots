#' Countsfile
#'
#' Count gene expression data for testing functions.
#' @docType data
#' @description
#' Dataset use to test final functions.
#' @format ## `countsfile` Contains 20 variables and 24500 observations
#' \describe{
#' \item{variables}{MeCP2 WT and KO groups}
#' \item{observations}{gene count levels}
#' }
#' @return dataset
"countsfile"

#' degsfile
#'
#' Differential Gene Expression (DEGs) data for testing functions.
#' @docType data
#' @description
#' Dataset use to test final functions.
#' @format ## `degsfile` Contains 4 variables and 13717 observations
#' \describe{
#' \item{logFC}{log Fold Change levels}
#' \item{logCPM}{log counts per million}
#' \item{PValue}{probability under null hypothesis}
#' \item{FDR}{False Discovery Rate}
#' }
#' @return dataset
"degsfile"

#' refseq
#'
#' Reference sequence data for testing functions.
#' @docType data
#' @description
#' Dataset use to test final functions.
#' @format ## `refseq` Contains 6 variables and 35976 observations
#' \describe{
#' \item{chr}{chromosome number}
#' \item{tx.start}{start location}
#' \item{tx.end}{end location}
#' \item{strand}{positive or negative}
#' \item{gene.name}{gene name}
#' \item{gene.length}{gene length in nucleotides}
#' }
#' @return dataset
"refseq"

#' mCA
#'
#' mCA data for testing functions.
#' @docType data
#' @description
#' Dataset use to test final functions.
#' @format ## `mCA` Contains 6 variables and 17916 observations
#' \describe{
#' \item{gene.name}{gene name}
#' \item{mCA}{methylated Cytosine-Adenine}
#' \item{CA}{Cytosine-Adenine}
#' \item{mCG}{methylated Cytosine-Guanine}
#' \item{CG}{Cytosine-Guanine}
#' \item{gene.length}{gene length in nucleotides}
#' }
#' @return dataset
"mCA"
