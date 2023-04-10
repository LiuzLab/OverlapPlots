#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 27th Nov 2019
#
# Program is used for:
# 1. Long and Short Gene analysis for Lisa Boxer Dataset
#############################################################################

## packages 
rm(list = ls())
setwd("~/Desktop/LongGene_Part-II/boxer/src/")
source("libraries.R")
library(ggpmisc)

## load the workspace
load(file = "../results/wholeCell.mecp2.Rdata")
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")
dat.annot <- data.frame(wholeCell.mecp2$counts) %>% rownames_to_column(var = "gene.name")
dat.annot <- inner_join(x = dat.annot, y = refseq, by ="gene.name")
dim(dat.annot)

## calculate the logFC between KO and WT
log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, 
                                           A.samples = c(2:11), 
                                           B.samples = c(12:21))

## logFC from boxer's DEGs list 
degs.dat <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_RNA-seq.txt",
                       sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
degs.dat <- degs.dat[degs.dat$FDR < 0.05,]
dim(degs.dat)
log2FC.comp <- inner_join(log2FC.length.KO[,c(1:21,30,26)],
                          degs.dat %>% rownames_to_column(var = "gene.name"), 
                          by = "gene.name")
eqs <- y ~ x
p1 <-  ggplot(data = log2FC.comp, aes(x = logFC.crude, y = logFC)) + 
        geom_point(size = 2) + geom_jitter() + theme_bw() +
        xlab("Log2FC -- DESeq2") + ylab("Log2FC -- edgeR") + 
        #geom_abline(col = "red", linetype = "dashed") +
        geom_smooth(method='lm', formula = eqs, col = "red", linetype = "dashed") +
        stat_poly_eq(formula = eqs, 
                    aes(label = paste(..eq.label..," ",..rr.label.., 
                                      sep = "~~~"), size = 20), parse = TRUE) +
        theme(plot.title = element_text(size = 18, face = "bold"),
              axis.title.y = element_text(size = 24, colour = "black",face = "bold"),
              axis.title.x = element_text(size = 24, colour = "black",face = "bold"),
              axis.text.y = element_text(size = 22, colour = "black",face = "bold"),
              axis.text.x = element_text(size = 22, colour = "black",face = "bold"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
ggsave(filename = "../results/Figures_Ver1/Supp1A.png", plot = p1, 
       width = 7, height = 6, dpi = 300, type = "cairo")

lm_fit <- lm(logFC.crude ~ logFC, data = log2FC.comp)
summary(lm_fit)
summary(lm_fit)$r.squared
cor(log2FC.comp$logFC, log2FC.comp$logFC.crude)

## Scatter Plot (Boxer vs ours analysis)
dat1 <- log2FC.comp[,c("gene.name","logFC", "FDR", "gene.length")]
Plot.Scatter(dat = dat1, log2FC = log2(1), comp.between = "Boxer's Whole Cell KO/WT")

dat2 <- inner_join(x = wholeCell.mecp2$results[,c(1,3,7)], 
                   y = log2FC.comp[,c("gene.name","gene.length")], 
                   by = c("gene" ="gene.name"))
Plot.Scatter(dat = dat2, log2FC = log2(1.20), pval = 0.05, comp.between = "Whole Cell KO/WT")

## data frame required for long gene plot
dat.annot <- log2FC.comp[,c(1:21,24,23)]
head(log2FC.comp)
colnames(log2FC.comp)[22] <- "logFC.crude"

## all the combinations now:
num.control <- 10
ind.WT <- c(1:10)
ind.KO <- c(11:20)
control.vec <- c(2:(num.control+1))
tot.comb <- combn(control.vec, 5)
iter <- ncol(tot.comb)/2
dir1 <- "../results/DEGs/whole_cell/"

idx1 <- c(2,3,4,7,11)
idx2 <- c(5,6,8,9,10)

for(i in 1:8){
  #idx1 <- c(2,3,4,6,7) #tot.comb[,i] 
  #idx2 <- c(5,8,9,10,11)  #control.vec[which(!control.vec %in% idx1)]
  cat("Combination",i,"include control groups as")
  cat("\n\tGroup 1:",colnames(log2FC.length.WT)[idx1],"and","\n\tGroup 2:",
      colnames(log2FC.length.WT)[idx2],"\n\n")

  ## logFC calculation
  log2FC.length.WT <- dat.annot[,c(1:(num.control+1),23)]
  log2FC.length.WT$comp.mat <- apply(dat.annot[,c(idx1,idx2)], 1,
                                     function(r){log2((mean(r[1:5])+1)/(mean(r[6:10])+1))})
  log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat")],
                              y = log2FC.length.KO[,c("gene.name","logFC.crude",
                                                      "gene.length")],
                              by = "gene.name")
  p3 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)],
                            comp.between1 = "(WT/WT)",
                            comp.between2 = "(KO/WT)", bin.size = 60, shift.size = 6)
  p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.25,0.25)),
                  p3$plot2 + coord_cartesian(ylim = c(0,30)),ncol = 1, align = 'v') +
                  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  filename <- paste0(dir1,"plot",i,".png")
  plot(p3)
  ggsave(plot = p3, filename = filename, width = 10, height = 8,
         dpi = 300, type = "cairo")
}

## problems or to do things:
## 1. specify the WT and KO sample indexes (done)
## 2. Add number of long bins that are significant and
## total number of long bins and genes covered; 
## shift size is very small (done)
## 3. create separate files 

