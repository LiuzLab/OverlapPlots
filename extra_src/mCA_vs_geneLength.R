#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 04th June 2020
#
# Program is used for:
# 1. mCA vs gene length plot
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/boxer/src/")
source("libraries.R")
source("overlap_wrappers.R")

## functions
lm_eqn <- function(df){
  m <- lm(mCA.CA ~ gene.length, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  #return(summary(m))
  as.character(as.expression(eq));
}

## load the workspace
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")
mCA <- data.frame(readRDS("../dat/CAperGene_10wk_boxer.RDS"), stringsAsFactors = F)
mCA <- mCA[mCA$CA >= 5,]
mCA <- mCA[mCA$gene.length >= 4500,]
mCA$mCA.CA <- mCA$mCA/mCA$CA
mCA1 <- mCA[,c(1,6,7)]
mCA1$gene.length <- mCA1$gene.length/1e3
mCA1$gl.type <- ifelse(mCA1$gene.length >= 100, "Long Gene", "Short Gene")

## scatter plot -- mCA vs gene length
formule <- y ~ x
ggplot(data = mCA1, aes(x = gene.length, y = mCA.CA, color = gl.type)) + 
  scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
  geom_point(alpha = 0.5, size = 0.25) + ylim(0,0.1) + 
  ylab("Gene body mCA/CA") + xlab("Gene Length in Kb ") +
  geom_smooth(method = lm, formula = formule, color = "black", size = 1, 
              linetype = "dashed") + facet_wrap(~ gl.type, dir="v") + theme_bw() + 
  theme(axis.title.y = element_text(size = 24, face = "bold", color = "black"),
        axis.text.y = element_text(size = 24, face = "bold", color = "black"),
        axis.title.x = element_text(size = 24, face = "bold", color = "black"),
        axis.text.x = element_text(size = 24, face = "bold", color = "black"),
        strip.text.x = element_text(size = 24, face = "bold", colour = "red"),
        legend.position="none")

summary(lm(mCA.CA ~ gene.length, 
           mCA1[mCA1$gl.type %in% "Short Gene",c(2,3)]))
lm_eqn(mCA1[mCA1$gl.type %in% "Short Gene",c(2,3)])

summary(lm(mCA.CA ~ gene.length,
           mCA1[mCA1$gl.type %in% "Long Gene",c(2,3)]))
lm_eqn(mCA1[mCA1$gl.type %in% "Long Gene",c(2,3)])

summary(lm(mCA.CA ~ gene.length, mCA1[,c(2,3)]))
lm_eqn(mCA1[,c(2,3)])

## DEGs

####################################
#
## whole cell dataset
#
####################################

## mCA analysis for KO/WT
mCA1 <- mCA[,c(1,6,7)]
load(file = "../results/wholeCell.KO.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(wholeCell.KO$counts[,c(1:10)])
res1 <- overlap_degs_mCA_wrapper(degs.dat, wholeCell.KO$counts[,c(1:10)],
                                 refseq = mCA1, WT1.idx = c(1,2,3,6,10), 
                                 WT2.idx = c(4,5,7,8,9), bin.size = 60,
                                 shift.size = 6, methyl.type = "mCA/CA", 
                                 degs = T)
p1 <- (res1$plot1 + coord_cartesian(ylim = c(-0.25,0.25))) / (res1$plot2)
rm(wholeCell.KO, degs.dat)

## mCA analysis for R306C/WT
load(file = "../results/wholeCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(wholeCell.R306C$counts[,c(1:10)])
res2 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = wholeCell.R306C$counts[,c(1:10)],
                                 refseq = mCA1, WT1.idx = c(7,8,9,10),
                                 WT2.idx = c(1,2,3,4,5,6), bin.size = 40, 
                                 shift.size = 4, methyl.type = "mCA/CA")
p2 <- (res2$plot1 + coord_cartesian(ylim = c(-0.25,0.25))) / (res2$plot2)
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
res3 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = nuclearCell.KO$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1,2,3), WT2.idx = c(4:10),
                                 bin.size = 40, shift.size = 4, 
                                 methyl.type = "mCA/CA")
p3 <- (res3$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res3$plot2)
rm(nuclearCell.KO, degs.dat)

## overlap plot R306C/WT
load(file = "../results/nuclearCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_nuclear_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(nuclearCell.R306C$counts[,c(1:10)])
res4 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = nuclearCell.R306C$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1,2,7,8,9), WT2.idx = c(3,4,5,6,10), 
                                 bin.size = 40, shift.size = 4,
                                 methyl.type = "mCA/CA")
p4 <- (res4$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res4$plot2)
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
res5 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = chromatinCell.KO$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1:4), WT2.idx = c(5:10), 
                                 bin.size = 40, shift.size = 4, 
                                 methyl.type = "mCA/CA")
p5 <- (res5$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res5$plot2)
rm(chromatinCell.KO, degs.dat)

## overlap plot R306C/WT whole cell
load(file = "../results/chromatinCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_chromatin_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(chromatinCell.R306C$counts[,c(1:10)])
res6 <- overlap_degs_mCA_wrapper(degs.dat, chromatinCell.R306C$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1,2,3,4,7,8,9), WT2.idx = c(5,6,10),
                                 bin.size = 40, shift.size = 4, 
                                 methyl.type = "mCA/CA")
p6 <- (res6$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res6$plot2)
rm(chromatinCell.R306C, degs.dat)

## combined mCA plot
(p1 | p3 | p5)/(p2 | p4 | p6)
((res1$plot2 + coord_cartesian(ylim = c(0,0.08), xlim = c(0,750))) | 
  (res3$plot5 + coord_cartesian(ylim = c(0,0.08), xlim = c(0,750))) | 
  (res5$plot5 + coord_cartesian(ylim = c(0,0.08), xlim = c(0,750)))) / 
  ((res2$plot5 + coord_cartesian(ylim = c(0,0.08), xlim = c(0,750))) | 
  (res4$plot5 + coord_cartesian(ylim = c(0,0.08), xlim = c(0,750))) | 
  (res6$plot5 + coord_cartesian(ylim = c(0,0.08), xlim = c(0,750))))

## all the genes

####################################
#
## whole cell dataset
#
####################################

## mCA analysis for KO/WT
mCA1 <- mCA[,c(1,6,7)]
load(file = "../results/wholeCell.KO.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(wholeCell.KO$counts[,c(1:10)])
res1 <- overlap_degs_mCA_wrapper(degs.dat, wholeCell.KO$counts[,c(1:10)],
                                 refseq = mCA1, WT1.idx = c(1,2,3,6,10), 
                                 WT2.idx = c(4,5,7,8,9), bin.size = 200,
                                 shift.size = 40, methyl.type = "mCA/CA", 
                                 degs = F)
p1 <- (res1$plot1 + coord_cartesian(ylim = c(-0.25,0.25))) / (res1$plot2)
rm(wholeCell.KO, degs.dat)

## mCA analysis for R306C/WT
load(file = "../results/wholeCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(wholeCell.R306C$counts[,c(1:10)])
res2 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = wholeCell.R306C$counts[,c(1:10)],
                                 refseq = mCA1, WT1.idx = c(7,8,9,10),
                                 WT2.idx = c(1,2,3,4,5,6), bin.size = 200, 
                                 shift.size = 40, methyl.type = "mCA/CA", 
                                 degs = F)
p2 <- (res2$plot1 + coord_cartesian(ylim = c(-0.25,0.25))) / (res2$plot2)
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
res3 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = nuclearCell.KO$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1,2,3), WT2.idx = c(4:10),
                                 bin.size = 200, shift.size = 40, 
                                 methyl.type = "mCA/CA", degs = F)
p3 <- (res3$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res3$plot2)
rm(nuclearCell.KO, degs.dat)

## overlap plot R306C/WT
load(file = "../results/nuclearCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_nuclear_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(nuclearCell.R306C$counts[,c(1:10)])
res4 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = nuclearCell.R306C$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1,2,7,8,9), WT2.idx = c(3,4,5,6,10), 
                                 bin.size = 200, shift.size = 40,
                                 methyl.type = "mCA/CA", degs = F)
p4 <- (res4$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res4$plot2)
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
res5 <- overlap_degs_mCA_wrapper(degs.dat, count.dat = chromatinCell.KO$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1:4), WT2.idx = c(5:10), 
                                 bin.size = 200, shift.size = 40, 
                                 methyl.type = "mCA/CA", degs = F)
p5 <- (res5$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res5$plot2)
rm(chromatinCell.KO, degs.dat)

## overlap plot R306C/WT whole cell
load(file = "../results/chromatinCell.R306C.Rdata")
degs.dat <- read.table("../dat/publication/DEGs/R306C-WT_chromatin_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(chromatinCell.R306C$counts[,c(1:10)])
res6 <- overlap_degs_mCA_wrapper(degs.dat, chromatinCell.R306C$counts[,c(1:10)], 
                                 mCA1, WT1.idx = c(1,2,3,4,7,8,9), WT2.idx = c(5,6,10),
                                 bin.size = 200, shift.size = 40, 
                                 methyl.type = "mCA/CA", degs = F)
p6 <- (res6$plot1 + coord_cartesian(ylim = c(-0.4,0.4))) / (res6$plot2)
rm(chromatinCell.R306C, degs.dat)

## combined mCA plot
(p1 | p3 | p5)/(p2 | p4 | p6)
((res1$plot5 + coord_cartesian(ylim = c(0.02,0.05), xlim = c(0,750))) | 
    (res3$plot5 + coord_cartesian(ylim = c(0.02,0.05), xlim = c(0,750))) | 
    (res5$plot5 + coord_cartesian(ylim = c(0.02,0.05), xlim = c(0,750)))) / 
  ((res2$plot5 + coord_cartesian(ylim = c(0.02,0.05), xlim = c(0,750))) | 
     (res4$plot5 + coord_cartesian(ylim = c(0.02,0.05), xlim = c(0,750))) | 
     (res6$plot5 + coord_cartesian(ylim = c(0.02,0.05), xlim = c(0,750))))

