#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 18th May 2020
#
# Program is used for:
# 1. Long and Short Gene analysis for all the Lisa Boxer Dataset
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/boxer/src/")
source("libraries.R")
source("overlap_wrappers.R")

## load the workspace
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")

####################################
#
## whole cell dataset
#
####################################

## overlap plot KO/WT
load(file = "../results/wholeCell.KO.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(wholeCell.KO$counts[,c(1:10)])
res1 <- overlap_degs_wrapper(degs.dat, count.dat = wholeCell.KO$counts[,c(1:10)], 
                             refseq, WT1.idx = c(1,2,3,6,10), WT2.idx = c(4,5,7,8,9),
                             bin.size = 60, shift.size = 6)
res1$plot
wholeCell.KO$plots$PCAplot_highvar
rm(wholeCell.KO, degs.dat)

## overlap plot R306C/WT
load(file = "../results/wholeCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(wholeCell.R306C$counts[,c(1:10)])
res2 <- overlap_degs_wrapper(degs.dat, count.dat = wholeCell.R306C$counts[,c(1:10)], 
                             refseq, WT2.idx = c(1,2,3,4,5,6), WT1.idx = c(7,8,9,10),
                             bin.size = 40, shift.size = 4)
wholeCell.R306C$plots$PCAplot_highvar
res2$plot
rm(wholeCell.R306C, degs.dat)

####################################
#
## nuclear dataset
#
####################################

## overlap plot KO/WT
load(file = "../results/nuclearCell.KO.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/KO-WT_nuclear_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(nuclearCell.KO$counts[,c(1:10)])
res3 <- overlap_degs_wrapper(degs.dat, count.dat = nuclearCell.KO$counts[,c(1:10)], 
                             refseq, WT1.idx = c(1:3), WT2.idx = c(4:10),
                             bin.size = 40, shift.size = 4)
nuclearCell.KO$plots$PCAplot_highvar
res3$plot
rm(nuclearCell.KO, degs.dat)

## overlap plot R306C/WT
load(file = "../results/nuclearCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_nuclear_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(nuclearCell.R306C$counts[,c(1:10)])
res4 <- overlap_degs_wrapper(degs.dat, count.dat = nuclearCell.R306C$counts[,c(1:10)], 
                             refseq, WT2.idx = c(1,2,7,8,9), WT1.idx = c(3,4,5,6,10), 
                             bin.size = 40, shift.size = 4)
nuclearCell.R306C$plots$PCAplot_highvar
res4$plot
rm(nuclearCell.R306C, degs.dat)

####################################
#
## chromatin associated dataset
#
####################################

## overlap plot KO/WT whole cell
load(file = "../results/chromatinCell.KO.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/KO-WT_chromatin_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(chromatinCell.KO$counts[,c(1:10)])
res5 <- overlap_degs_wrapper(degs.dat, count.dat = chromatinCell.KO$counts[,c(1:10)], 
                             refseq, WT1.idx = c(1:4), WT2.idx = c(5:10), 
                             bin.size = 40, shift.size = 4)
res5$plot
chromatinCell.KO$plots$PCAplot_highvar
rm(chromatinCell.KO, degs.dat)

## overlap plot R306C/WT whole cell
load(file = "../results/chromatinCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_chromatin_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(chromatinCell.R306C$counts[,c(1:10)])
res6 <- overlap_degs_wrapper(degs.dat, count.dat = chromatinCell.R306C$counts[,c(1:10)], 
                             refseq, WT2.idx = c(1,2,3,4,7,8,9), WT1.idx = c(5,6,10),
                             bin.size = 40, shift.size = 4)
res6$plot
chromatinCell.R306C$plots$PCAplot_highvar
rm(chromatinCell.R306C, degs.dat)

## saving the plots
res1$plot
ggsave(plot = res1$plot, 
       filename = "../results/Figures_Ver1/Fig1B_whole-cell_KO-WT_DEGs.png",
       width = 8, height = 6, dpi = 300, type = "cairo")

res2$plot
ggsave(plot = res2$plot, 
       filename = "../results/Figures_Ver1/Fig1D_whole-cell_R306C-WT_DEGs.png",
       width = 8, height = 6, dpi = 300, type = "cairo")

res3$plot
ggsave(plot = res3$plot, 
       filename = "../results/Figures_Ver1/Supp2B_nuclear-cell_KO-WT_DEGs.png",
       width = 8, height = 6, dpi = 300, type = "cairo")

res4$plot
ggsave(plot = res5$plot,
       filename = "../results/Figures_Ver1/Supp3B_chromatin_KO-WT_DEGs.png",
       width = 8, height = 6, dpi = 300, type = "cairo")

res5$plot
ggsave(plot = res5$plot,
       filename = "../results/Figures_Ver1/Supp3B_chromatin_KO-WT_DEGs.png",
       width = 8, height = 6, dpi = 300, type = "cairo")


q1 <- (res1$plot | res3$plot | res5$plot)/(res2$plot | res4$plot | res6$plot)
ggsave(plot = q1, 
       filename = "../results/Figures_Ver1/Boxer_LongGeneTrendPlots.tiff",
       width = 24, height = 12, dpi = 300, type = "cairo")

