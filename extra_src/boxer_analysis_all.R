#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 27th Dec 2019
#
# Program is used for:
# 1. Long and Short Gene analysis for all the Lisa Boxer Dataset
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/boxer/src/")
source("libraries.R")

## annotation dataset
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))

## KO/WT whole cell dataset --> 2922 DEGs
dat <- read.table("../dat/publication/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
wholeCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
ggsave(plot = wholeCell.KO$plot$PCAplot_highvar, 
       filename = "../results/Figures_Ver1/whole-cell_KO-WT_PCA.png",
       width = 8, height = 6, dpi = 300, type = "cairo")
#save(wholeCell.KO, file = "../results/wholeCell.KO.Rdata")
rm(dat)

## KO/WT nuclear dataset --> 2124 DEGs
dat <- read.table("../dat/publication/GSE128178_10WT_10MeCP2_KO_nuclear_RNAseq_gene_body_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
nuclearCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
#save(nuclearCell.KO, file = "../results/nuclearCell.KO.Rdata")
rm(dat)

## KO/WT chromatin dataset --> 1797 DEGs
dat <- read.table("../dat/publication/GSE128178_10WT_10MeCP2_KO_chromatin_associated_RNAseq_gene_body_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
chromatinCell.KO <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
#save(chromatinCell.KO, file = "../results/chromatinCell.KO.Rdata")
rm(dat)

## R306C/WT whole cell dataset --> 1271 DEGs
genotypes <- factor(c(rep("WT", 10), rep("MUT", 10)), levels = c("MUT", "WT"))
dat <- read.table("../dat/publication/GSE128178_10WT_10R306C_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
wholeCell.R306C <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
#save(wholeCell.R306C, file = "../results/wholeCell.R306C.Rdata")
rm(dat)

## R306C/WT nuclear dataset --> 1010 DEGs 
dat <- read.table("../dat/publication/GSE128178_10WT_10R306C_nuclear_RNAseq_gene_body_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
nuclearCell.R306C <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
#save(nuclearCell.R306C, file = "../results/nuclearCell.R306C.Rdata")
rm(dat)

## R306C/WT chromatin dataset --> 842 DEGs  
dat <- read.table("../dat/publication/GSE128178_10WT_10R306C_chromatin_associated_RNAseq_gene_body_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
chromatinCell.R306C <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
#save(chromatinCell.R306C, file = "../results/chromatinCell.R306C.Rdata")
rm(dat)

####################################
#
## whole cell dataset
#
####################################

plot_grid(wholeCell.KO$plots$PCAplot_highvar, 
          wholeCell.R306C$plots$PCAplot_highvar, ncol = 1, align = "v")

## overlap plot
source("overlap_wrappers.R")
res1 <- overlap_wrapper(dat = wholeCell.KO$counts, refseq = refseq, 
                       KO.idx = c(11:20), WT.idx = c(1:10), 
                       WT1.idx = c(1,2,3,6,10), WT2.idx = c(4,5,7,8,9))
res2 <- overlap_wrapper(dat = wholeCell.R306C$counts, refseq = refseq, 
                        KO.idx = c(11:20), WT.idx = c(1:10), 
                        WT1.idx = c(2,3,4,5,6), WT2.idx = c(1,7,8,9,10))
res1$plot
res2$plot

####################################
#
## nuclear dataset
#
####################################

plot_grid(nuclearCell.KO$plots$PCAplot_highvar, 
          nuclearCell.R306C$plots$PCAplot_highvar, ncol = 1, align = "v")

## overlap plot
res3 <- overlap_wrapper(dat = nuclearCell.KO$counts, refseq = refseq, 
                        KO.idx = c(11:20), WT.idx = c(1:10), 
                        WT1.idx = c(1,2,3,7,9), WT2.idx = c(4,5,6,8,10))
res4 <- overlap_wrapper(dat = nuclearCell.R306C$counts, refseq = refseq, 
                        KO.idx = c(11:20), WT.idx = c(1:10), 
                        WT1.idx = c(1,2,7,9,10), WT2.idx = c(3,4,5,6,8))
res3$plot
res4$plot


####################################
#
## chromatin associated dataset
#
####################################

plot_grid(chromatinCell.KO$plots$PCAplot_highvar, 
          chromatinCell.R306C$plots$PCAplot_highvar, ncol = 1, align = "v")

## overlap plot
res5 <- overlap_wrapper(dat = chromatinCell.KO$counts, refseq = refseq, 
                        KO.idx = c(11:20), WT.idx = c(1:10), 
                        WT1.idx = c(1:5), WT2.idx = c(6:10))
res6 <- overlap_wrapper(dat = chromatinCell.R306C$counts, refseq = refseq, 
                        KO.idx = c(11:20), WT.idx = c(1:10), 
                        WT1.idx = c(1,3,4,8,9), WT2.idx = c(2,5,6,7,10))
res5$plot
res6$plot




# dat.annot <- data.frame(wholeCell.KO$counts,
#                         gene.name = rownames(wholeCell.KO$counts),
#                         stringsAsFactors = F)
# dat.annot <- inner_join(x = dat.annot, y = refseq, by = "gene.name")
# 
# 
# log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, 
#                                            A.samples = c(1:10), 
#                                            B.samples = c(11:20))
# 
# ## random sampling into 2 groups based on WT samples
# cols.samples <- c(1:10)
# idx1 <- c(1,2,3,6,10)
# idx2 <- c(4,5,7,8,9)
# cat("Control Group 1:",idx1,"\n")
# cat("Control Group 2:",idx2,"\n")
# log2FC.length.WT <- dat.annot[,c(1:10,21)]
# log2FC.length.WT$comp.mat <- apply(dat.annot[,c(idx1,idx2)], 1,
#                                    function(r){log2((mean(r[1:5])+1)/(mean(r[6:10])+1))})
# log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat")],
#                             y = log2FC.length.KO[,c("gene.name","logFC.crude",
#                                                     "gene.length")], by = "gene.name")
# p1 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)],
#                           comp.between1 = "(WT/WT)",
#                           comp.between2 = "(KO/WT)")
# p1 <- plot_grid(p1$plot1 + coord_cartesian(ylim = c(-0.25,0.25)),
#                 p1$plot2 + coord_cartesian(ylim = c(0,50)),ncol = 1, align = 'v') +
#   theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
# p1
# 
