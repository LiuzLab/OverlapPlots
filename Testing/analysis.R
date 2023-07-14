#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Updated Date: 5th April 2023
#
# Program is used for:
# 1. Overlap analyses for all the Boxer et al Dataset
#############################################################################

## packages
rm(list = ls())
setwd("C:/Users/danie/Documents/longene2/R")
#

## annotation dataset
load(file = "../dat-info/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))

##################################
#
## KO/WT whole cell dataset
#
##################################

dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
mat <- wholeCell.KO$results[,c(8:27,1,3)]
colnames(mat)[21] <- "gene.name"

## same cluster sizes
grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
res1 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 200,
                        shift.size = 40)
res1$plot

### Using the DEGs list
degs.dat <- read.table("../dat/DEGs/KO-WT_whole-cell_RNA-seq.txt", sep = "\t",
                       stringsAsFactors=F, header = TRUE, row.names = 1)
degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
                       refseq, by = "gene.name")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1), comp.between = "")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1.2), comp.between = "")
degs.dat <- inner_join(x = wholeCell.KO$results[wholeCell.KO$results$gene %in%
                                                degs.dat$gene.name,c(8:27,1,3)],
                       y = degs.dat[,c("gene.name", "logFC")],
                       by = c("gene" = "gene.name"))
mat <- degs.dat[,c(1:21,23)]
colnames(mat)[21] <- "gene.name"

### logFC from Degs from edgeR (Boxer et al.) and
### much closer to the Boxer et al. paper
colnames(mat)[22] <- "log2FoldChange"
res2 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 60,
                        shift.size = 6, shrink_lfc = T)
res2$plot

### logFC from Degs from DESeq2
mat <- degs.dat[,c(1:22)]
colnames(mat)[21] <- "gene.name"
res3 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 60,
                        shift.size = 6, shrink_lfc = T)
res3$plot
rm(dat, wholeCell.KO, mat, grp.idx, degs.dat)

##################################
#
## KO/WT nuclear dataset
#
##################################

dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_nuclear_RNAseq_gene_body_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
nuclear.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
mat <- nuclear.KO$results[,c(8:27,1,3)]
colnames(mat)[21] <- "gene.name"
grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
    message("Cluster size is not equal, therefore run same size k-means variation!")
    grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
}
res4 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 200,
                        shift.size = 40)
res4$plot

### Using the DEGs list
degs.dat <- read.table("../dat/DEGs/KO-WT_nuclear_RNA-seq.txt", sep = "\t",
                       stringsAsFactors=F, header = TRUE, row.names = 1)
degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
                       refseq, by = "gene.name")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1), comp.between = "")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1.2), comp.between = "")
degs.dat <- inner_join(x = nuclear.KO$results[nuclear.KO$results$gene %in%
                                                    degs.dat$gene.name,c(8:27,1,3)],
                       y = degs.dat[,c("gene.name", "logFC")],
                       by = c("gene" = "gene.name"))
mat <- degs.dat[,c(1:21,23)]
colnames(mat)[21] <- "gene.name"

### logFC from Degs from edgeR (Boxer et al.) and
### much closer to the Boxer et al. paper
colnames(mat)[22] <- "log2FoldChange"
res5 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 40,
                        shift.size = 4, shrink_lfc = T)
res5$plot

### logFC from Degs from DESeq2
mat <- degs.dat[,c(1:22)]
colnames(mat)[21] <- "gene.name"
res6 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 40,
                        shift.size = 4, shrink_lfc = T)
res6$plot
rm(dat, nuclear.KO, mat, grp.idx, degs.dat)

#######################################
#
## KO/WT chromatin dataset
#
#######################################

dat <- read.table("../dat/counts/GSE128178_10WT_10MeCP2_KO_chromatin_associated_RNAseq_gene_body_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
chromatin.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
mat <- chromatin.KO$results[,c(8:27,1,3)]
colnames(mat)[21] <- "gene.name"
grp.idx <- WTgrp_kmeans(control_mat = mat[,1:10])
if(length(grp.idx$WT.idx1) != length(grp.idx$WT.idx2)){
    message("Cluster size is not equal, therefore run same size k-means variation!")
    grp.idx <- WTgrp_kmeans_eqSize(control_mat = mat[,1:10])
}
res7 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 200,
                        shift.size = 40)
res7$plot

### Using the DEGs list
degs.dat <- read.table("../dat/DEGs/KO-WT_chromatin_RNA-seq.txt", sep = "\t",
                       stringsAsFactors=F, header = TRUE, row.names = 1)
degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
cat("Number of DEGs =",dim(degs.dat)[1],"\n\n")
degs.dat <- inner_join(degs.dat %>% rownames_to_column(var = "gene.name"),
                       refseq, by = "gene.name")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1), comp.between = "")
Plot.Scatter(dat = degs.dat[,c("gene.name","logFC", "FDR", "gene.length")],
             log2FC = log2(1.2), comp.between = "")
degs.dat <- inner_join(x = chromatin.KO$results[chromatin.KO$results$gene %in%
                                                  degs.dat$gene.name,c(8:27,1,3)],
                       y = degs.dat[,c("gene.name", "logFC")],
                       by = c("gene" = "gene.name"))
mat <- degs.dat[,c(1:21,23)]
colnames(mat)[21] <- "gene.name"

### logFC from Degs from edgeR (Boxer et al.) and
### much closer to the Boxer et al. paper
colnames(mat)[22] <- "log2FoldChange"
res8 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 40,
                        shift.size = 4, shrink_lfc = T)
res8$plot

### logFC from Degs from DESeq2
mat <- degs.dat[,c(1:22)]
colnames(mat)[21] <- "gene.name"
res9 <- overlap_wrapper(dat = mat, refseq = refseq, KO.idx = c(11:20),
                        WT.idx = c(1:10), WT1.idx = grp.idx$WT.idx1,
                        WT2.idx = grp.idx$WT.idx2, bin.size = 40,
                        shift.size = 4, shrink_lfc = T)
res9$plot
rm(dat, chromatin.KO, mat, grp.idx, degs.dat)

### overlap_degs_wrapper
# res2 <- overlap_degs_wrapper(degs.dat, count.dat = wholeCell.KO$counts,
#                              refseq, WT1.idx = grp.idx$WT.idx1,
#                              WT2.idx = grp.idx$WT.idx2, bin.size = 60,
#                              shift.size = 6)
# res2$plot
#
