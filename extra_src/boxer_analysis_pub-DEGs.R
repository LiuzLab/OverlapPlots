#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 17th Dec 2019
#
# Program is used for:
# 1. Lisa Boxer DEGs Dataset vs DESeq analysis
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/boxer/src/")
source("libraries.R")

## datasets
load(file = "../results/wholeCell.mecp2.Rdata")
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")

## degs from boxer
boxer.degs <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_RNA-seq.txt",
                         sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
boxer.degs <- boxer.degs %>% tibble::rownames_to_column(var = "gene.name")
boxer.degs <- inner_join(x = boxer.degs[,c(1,2,5)], 
                         y = refseq[,c("gene.name","gene.length")], 
                         by = "gene.name")
Plot.Scatter(dat = boxer.degs, log2FC = log2(1), pval = 0.05, 
             comp.between = "Whole Cell KO/WT (Lisa's analysis)")

## DESeq Degs
DEseq.degs <- inner_join(x = wholeCell.mecp2$results[,c(1,3,7)], 
                         y = refseq[,c("gene.name","gene.length")], 
                         by = c("gene"="gene.name"))
Plot.Scatter(dat = DEseq.degs, log2FC = log2(1), pval = 0.05, 
             comp.between = "Whole Cell KO/WT (DEseq)")

## Common DEGs and its distribution
DEseq.degs$Sig1 <- ifelse(DEseq.degs$log2FoldChange > 0 & DEseq.degs$padj < 0.05, "Up",
                          ifelse(DEseq.degs$log2FoldChange < 0 & DEseq.degs$padj < 0.05, "Down",
                                 "Not Sigficant"))
boxer.degs$Sig2 <- ifelse(boxer.degs$logFC > 0 & boxer.degs$FDR < 0.05, "Up",
                         ifelse(boxer.degs$logFC < 0 & boxer.degs$FDR < 0.05, "Down",
                                "Not Sigficant"))
table(DEseq.degs$Sig1)
table(boxer.degs$Sig2)

both.degs <- inner_join(x = DEseq.degs, y = boxer.degs, by = c("gene"="gene.name"))
sum(both.degs$Sig1 == "Up" & both.degs$Sig2 == "Up") 
sum(both.degs$Sig1 == "Down" & both.degs$Sig2 == "Down") 

sum(both.degs$Sig1 == "Up" & both.degs$Sig2 == "Not Sigficant") 
sum(both.degs$Sig1 == "Down" & both.degs$Sig2 == "Not Sigficant") 

sum(both.degs$Sig1 == "Not Sigficant" & both.degs$Sig2 == "Up") 
sum(both.degs$Sig1 == "Not Sigficant" & both.degs$Sig2 == "Down") 

## Degs analysis based on the PCA
## Group1 -- c(1,2,3,7); 
## Group2 -- c(4,5,6,9,10) 
source(file = "DESeqCalculation_WTcomp.R")

wholeCell.mecp2$plots$PCAplot_rld
idx1 <- c(1,2,3,6,10)
idx2 <- c(4,5,7,8,9)

dat <- read.table("../dat/publication/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
genotypes <- factor(c(rep("WT1", 5), rep("WT2", 5)), levels = c("WT1", "WT2"))
wholeCell.WT.mecp2 <- DESeqCalculation_WTcomp(dat = dat[,c(idx1,idx2)], 
                                       genotypes = genotypes, fc = 1)
wholeCell.WT.degs <- inner_join(x = wholeCell.WT.mecp2$results[,c(1,3,7)], 
                         y = refseq[,c("gene.name","gene.length")], 
                         by = c("gene"="gene.name"))

Plot.Scatter(dat = wholeCell.WT.degs, log2FC = log2(1), comp.between = "(WT1/WT2)")
Plot.Scatter(dat = wholeCell.WT.degs, log2FC = log2(1.1), comp.between = "(WT1/WT2)")
Plot.Scatter(dat = wholeCell.WT.degs, log2FC = log2(1.20), comp.between = "(WT1/WT2)")
Plot.Scatter(dat = wholeCell.WT.degs, log2FC = log2(1.30), comp.between = "(WT1/WT2)")


wholeCell.WT.mecp2$plots$sizefactorPlot + 
  ggrepel::geom_label_repel(aes(label = colData$samples), size = 6, fontface = "bold",
                           color = "black", box.padding = unit(0.35, "lines"),
                           point.padding = unit(0.2, "lines"))

wholeCell.mecp2$plots$sizefactorPlot + 
  ggrepel::geom_text_repel(aes(label = coData$samples), size = 6, fontface = "bold",
                           color = "black", box.padding = unit(0.35, "lines"),
                           point.padding = unit(0.2, "lines"))
save(wholeCell.WT.mecp2, file = "../results/WT_comp/wholeCell.WT.Rdata")




