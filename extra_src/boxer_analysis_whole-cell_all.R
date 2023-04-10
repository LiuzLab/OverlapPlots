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

## dataset -- 1
dat <- read.table("../dat/publication/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                            sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
head(dat)
genotypes <- factor(c(rep("WT", 10), rep("KO", 10)), levels = c("KO", "WT"))
wholeCell.mecp2 <- DESeqCalculation(dat = dat, genotypes = genotypes, fc = 1.15)
save(wholeCell.mecp2, file = "../results/wholeCell.mecp2.Rdata")

## Mecp2
ind <- which(rownames(wholeCell.mecp2$counts) == "Mecp2")  ## Mecp2, "NM_001110792"
plot.dat <- data.frame(labels = factor(colnames(dat[,c(11:20,1:10)]), 
                                       levels = colnames(dat[,c(11:20,1:10)])),
                       Normalized.Counts = as.vector(
                         as.matrix(log2(dat[ind,c(11:20,1:10)])+1)), 
                       genotypes = relevel(genotypes[c(11:20,1:10)],"KO"))
print(ggplot(plot.dat, aes(x = labels, y = Normalized.Counts, fill = genotypes)) +
        geom_bar(stat="identity", width=0.75, position = position_dodge(width=0.5)) +
        xlab("") + ylab("Normalized Counts") + ggtitle(paste("Barplot for Gene MeCP2")) + theme_bw() +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              legend.text = element_text(size = 12, face = "bold"),
              legend.title = element_text(size = 12, colour = "black",face = "bold"),
              axis.title.y= element_text(size = 14, colour = "black", face = "bold"),
              axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
              axis.text.x = element_blank(), axis.ticks.x = element_blank()))

## long gene annotation
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")
dat.annot <- data.frame(wholeCell.mecp2$counts) %>% rownames_to_column(var = "gene.name")
dat.annot <- inner_join(x = dat.annot, y = refseq, by ="gene.name")
dim(dat.annot)

## overlap plot
log2FC.length.KO <- logofMeans.between.A.B(dat = dat.annot, 
                                           A.samples = c(2:11), 
                                           B.samples = c(12:21))

## random sampling into 2 groups based on WT samples
cols.samples <- c(2:11)
idx1 <- c(2,3,4,7,11)
idx2 <- c(5,6,8,9,10)
cat("Control Group 1:",idx1,"\n")
cat("Control Group 2:",idx2,"\n")
log2FC.length.WT <- dat.annot[,c(1:11,26)]
log2FC.length.WT$comp.mat <- apply(dat.annot[,c(idx1,idx2)], 1,
                                   function(r){log2((mean(r[1:5])+1)/(mean(r[6:10])+1))})
log2FC.length <- inner_join(x = log2FC.length.WT[,c("gene.name","comp.mat")],
                            y = log2FC.length.KO[,c("gene.name","logFC.crude",
                                                    "gene.length")], by = "gene.name")
p3 <- overlay.gabels.plot(mat = log2FC.length[,c(2:4)],
                          comp.between1 = "(WT/WT)",
                          comp.between2 = "(KO/WT)")
p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.25,0.25)),
                p3$plot2 + coord_cartesian(ylim = c(0,50)),ncol = 1, align = 'v') +
                theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
p3

## all the combinations now:
num.control <- 10
ind.WT <- c(1:10)
ind.KO <- c(11:20)
control.vec <- c(2:(num.control+1))
tot.comb <- combn(control.vec, 5)
iter <- ncol(tot.comb)/2
dir1 <- "../results/All genes/whole_cell_check/"

for(i in 1:iter){
  idx1 <- tot.comb[,i]
  idx2 <- control.vec[which(!control.vec %in% idx1)]
  cat("\n\n")
  cat("Combination",i,"include control groups as ")
  cat("Group 1:",idx1-1,"and","Group 2:",idx2-1,"\n")

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
                            comp.between2 = "(KO/WT)", bin.size = 200, shift.size = 40)
  p3 <- plot_grid(p3$plot1 + coord_cartesian(ylim = c(-0.25,0.25)),
                  p3$plot2 + coord_cartesian(ylim = c(0,30)),ncol = 1, align = 'v') +
                  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
  filename <- paste0(dir1,"plot",i,".png")
  ggsave(plot = p3, filename = filename, width = 10, height = 8,
         dpi = 300, type = "cairo")
}
