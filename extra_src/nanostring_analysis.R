#############################################################################
# @Author: Ayush T. Raman
# Liu Lab, NRI, Baylor College of Medicine
# Date: 4th Dec 2019
#
# Program is used for:
# 1. Analysis of Boxer nanostring dataset
#############################################################################

rm(list = ls())

## functions and set the working directory
setwd("~/Desktop/LongGene_Part-II/boxer/src/")
source("libraries.R")

## Taking the normalized dataset from nSolver (as processed by Boxer et al.)
norm.table <- read.table(file = "../dat/publication/boxer_nanostring_KO-WT_normalized_nSolver.txt", 
                         header = F, sep = "\t", quote = "", fill = T,
                         stringsAsFactors = FALSE, na.strings = FALSE)
norm.table <- norm.table[-c(1,3,204:217),-c(2:7)]
col.names <- norm.table[1,c(1:21)]
col.names <- gsub(pattern = "(.*)_(WT|KO)([0-9]+)_(.*)",replacement = "\\2\\3", x = col.names)
col.names[1] <- "gene.name"
colnames(norm.table) <- col.names
norm.table <- norm.table[-1,]

## refseq file
load(file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")
rownames(refseq) <- refseq$gene.name
sum(norm.table$gene.name %in% refseq$gene.name)
norm.table[which(!norm.table$gene.name %in% refseq$gene.name),"gene.name"]

## change the inconsistent gene namess
idx <- which(norm.table$gene.name == "Epb4.1l4b")
norm.table$gene.name[idx] <- "Epb41l4b"
idx <- which(norm.table$gene.name == "Lphn2")
norm.table$gene.name[idx] <- "Adgrl2"
idx <- which(norm.table$gene.name == "Mar4_")
norm.table$gene.name[idx] <- "March4"
idx <- which(norm.table$gene.name == "MeCP2")
norm.table$gene.name[idx] <- "Mecp2"
idx <- which(norm.table$gene.name == "Ppap2a")
norm.table$gene.name[idx] <- "Plpp1"
idx <- which(norm.table$gene.name == "Sep6_")
norm.table$gene.name[idx] <- "Sept6"
idx <- which(norm.table$gene.name == "Zak (Map3k20)")
norm.table$gene.name[idx] <- "Zak"
sum(norm.table$gene.name %in% refseq$gene.name)

## PCA plot
dim(norm.table)
dat <- apply(norm.table[,c(2:21)], 2, function(r) as.numeric(as.character(r)))
rownames(dat) <- norm.table$gene.name
genotypes <- factor(rep(c("WT", "KO"), each = 10), levels = c("KO", "WT"))
p1 <- PCAplot(data = log2(dat+1), genotypes = genotypes)
p2 <- MDSplot(data = log2(dat+1), genotypes = genotypes)

## Comparison plot (RNA-seq and Nanostring)
degs.rnaseq <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_RNA-seq.txt",
                          sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(degs.rnaseq) <- paste0(colnames(degs.rnaseq),"_RNAseq")
degs.nano <- read.table("../dat/publication/DEGs/KO-WT_whole-cell_nanoString.txt",
                        sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
colnames(degs.nano)[c(5,6,7)] <- c("LogFC_nano","PValue_nano","FDR_nano")
degs.nano <- degs.nano[,c(1,3,5,6,7)]
degs.annot <- inner_join(degs.rnaseq %>% rownames_to_column(var = "gene.name"), 
                         degs.nano %>% rownames_to_column(var = "gene.name"), 
                         by = "gene.name")
degs.annot <- inner_join(x = degs.annot, y = refseq, by = "gene.name")
degs.annot$Sig <- ifelse(degs.annot$`FDR_RNAseq` < 0.05 & degs.annot$FDR_nano < 0.05, "Both",
                          ifelse(degs.annot$`FDR_RNAseq` < 0.05 & degs.annot$FDR_nano >= 0.05, "RNA-seq",
                          ifelse(degs.annot$`FDR_RNAseq`>= 0.05 & degs.annot$FDR_nano < 0.05, "Nanostring",
                                 "Not Significant")))
degs.annot$long.short <- ifelse(degs.annot$gene.length > 100e3, "Long", "Short")
degs.annot$length.sig <- paste0(degs.annot$long.short,"--",degs.annot$Sig)

p1 <- qplot(data = degs.annot,#[degs.annot$long.short == "Short",], 
              x = logFC_RNAseq, y = LogFC_nano, col = length.sig) + 
        geom_point(shape=1) + geom_jitter() + geom_abline(colour = "grey", size = 0.5) + 
        xlim(c(-1.5, 1)) + ylim(c(-1.5, 1)) + ggtitle("Long Genes") + theme_bw() +
        geom_hline(aes(yintercept=0), linetype="dotted") + 
        geom_vline(aes(xintercept=0), linetype="dotted") +        
        geom_hline(aes(yintercept=log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=log2(1.2)), linetype="dotted") +
        geom_hline(aes(yintercept=-log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=-log2(1.2)), linetype="dotted") +
        xlab(paste("LogFC RNA-Seq Values")) + ylab(paste("LogFC NanoString Values")) + 
        theme(legend.text = element_text(size = 12, face = "bold"), 
              legend.title = element_text(size = 12, colour = "black", face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
              plot.title = element_blank())

plot(degs.annot$logFC_RNAseq[degs.annot$long.short == "Long"])
abline(h = log2(1.2))
abline(h = log2(1/1.2))

plot(degs.annot$logFC_RNAseq[degs.annot$long.short == "Short"])
abline(h = log2(1.2))
abline(h = log2(1/1.2))

## Figure 3 
q1 <- qplot(y = degs.annot$logFC_RNAseq[degs.annot$long.short == "Long"],
            x = 1:length(degs.annot$logFC_RNAseq[degs.annot$long.short == "Long"])) +
        geom_point(size = 1) + theme_classic() + 
        xlab("Index") + ylab("Log2FC (Long Genes)") + 
        geom_hline(yintercept = log2(1.2), linetype = "dashed") + 
        geom_hline(yintercept = log2(1/1.2), linetype = "dashed") +
        theme(axis.title = element_text(size = 18, face = "bold"), legend.position="none", 
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              legend.text = element_text(size = 18, face = "bold"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 

q2 <- qplot(y = degs.annot$logFC_RNAseq[degs.annot$long.short == "Short"],
            x = 1:length(degs.annot$logFC_RNAseq[degs.annot$long.short == "Short"])) +
        geom_point(size = 1) + theme_classic() + 
        xlab("Index") + ylab("Log2FC (Short Genes)") + 
        geom_hline(yintercept = log2(1.2), linetype = "dashed") + 
        geom_hline(yintercept = log2(1/1.2), linetype = "dashed") + 
        theme(axis.title = element_text(size = 18, face = "bold"), legend.position="none", 
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              legend.text = element_text(size = 18, face = "bold"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
q1 | q2

## LogFC  diff plot
degs.annot$abs.logFC <- abs(degs.annot$logFC_RNAseq) - abs(degs.annot$LogFC_nano)
long.up <- sum(degs.annot$abs.logFC > 0 & degs.annot$long.short == "Long")
long.down <- sum(degs.annot$abs.logFC <=  0 & degs.annot$long.short == "Long")  
short.up <- sum(degs.annot$abs.logFC > 0 & degs.annot$long.short == "Short")
short.down <- sum(degs.annot$abs.logFC <=  0 & degs.annot$long.short == "Short")  

ggplot(data = degs.annot, aes(x = gene.length/1000, y = abs.logFC, col = long.short)) + 
  geom_point(size = 2) + geom_hline(aes(yintercept = 0), 
                                        colour="#FF0000", linetype="dashed", 
                                        size = 1) + 
  scale_x_continuous(trans = log10_trans(), breaks = c(1,10,100,1000)) +
  ylab(paste("abs(RNA-seq - Nanostring) LogFC (Classical)")) + 
  xlab(paste("Gene Length")) + theme_bw() + ylim(c(-0.5,0.5)) +
  annotate("text",x = 500,y = 0.35,label = long.up, size = 7, fontface = "bold") +
  annotate("text",x = 500,y = -0.35,label = long.down, size = 7, fontface = "bold") +
  annotate("text",x = 5,y = 0.35,label = short.up, size = 7, fontface = "bold") +
  annotate("text",x = 5,y = -0.35,label = short.down, size = 7, fontface = "bold") +
  theme(legend.position="none", 
        axis.title = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 24, face = "bold", color = "black"),
        axis.text.y = element_text(size = 24, face = "bold", color = "black"))

###################################
#
## Now compare shrunken logFC
#
###################################

load("../results/wholeCell.mecp2.Rdata")

## DEGs analysis
dds.nano <- DESeqDataSetFromMatrix(countData = as.matrix(round(dat)),
                                 colData = data.frame(row.names = colnames(dat),
                                                       genotypes = genotypes),
                                  design = ~ genotypes)
sizeFactors(dds.nano) <- rep(1, 20)
dds.nano <- DESeq(dds.nano, betaPrior = TRUE)
res.nano <- results(dds.nano, contrast = c("genotypes", "KO", "WT"))
res.seqnano <- inner_join(x = data.frame(wholeCell.mecp2$results), 
                          y = data.frame(res.nano) %>% rownames_to_column(var = "gene"),
                          by = "gene")
res.seqnano <-inner_join(x = res.seqnano,  y = refseq, by =c("gene"="gene.name"))
res.seqnano$abs.logFC <- res.seqnano$log2FoldChange.x - res.seqnano$log2FoldChange.y
res.seqnano$long.short <- ifelse(res.seqnano$gene.length > 100e3, "Long", "Short")

ggplot(data = res.seqnano, aes(x = gene.length/1000, y = abs.logFC, col = long.short)) + 
  geom_point(size = 2) + geom_hline(aes(yintercept = 0), 
                                    colour="#FF0000", linetype="dashed", 
                                    size = 1) + 
  scale_x_continuous(trans = log10_trans(), breaks = c(1,10,100,1000)) +
  ylab(paste("abs(RNA-seq - Nanostring) LogFC (Classical)")) + 
  xlab(paste("Gene Length")) + theme_bw() +
  theme(legend.position="none", 
        axis.title = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 24, face = "bold", color = "black"),
        axis.text.y = element_text(size = 24, face = "bold", color = "black"))

ggplot(res.seqnano, aes(x = log2FoldChange.x, y = res.seqnano$log2FoldChange.y)) + 
  geom_point(size = 2) + xlim(c(-2,2)) + ylim(c(-2,2))

cor(res.seqnano$log2FoldChange.x, res.seqnano$log2FoldChange.y)
cor(degs.annot$logFC_RNAseq, degs.annot$LogFC_nano)

res.seqnano$Sig <- ifelse(res.seqnano$padj.x < 0.05 & res.seqnano$padj.y < 0.05, "Both",
                         ifelse(res.seqnano$padj.x < 0.05 & res.seqnano$padj.y >= 0.05, "RNA-seq",
                                ifelse(res.seqnano$padj.x >= 0.05 & res.seqnano$padj.y < 0.05, "Nanostring",
                                       "Not Significant")))
res.seqnano$long.short <- ifelse(res.seqnano$gene.length > 100e3, "Long", "Short")
res.seqnano$length.sig <- paste0(res.seqnano$long.short,"--",degs.annot$Sig)

p2 <- qplot(data = res.seqnano[res.seqnano$long.short == "Long",], 
            x = log2FoldChange.x, y = log2FoldChange.y, col = length.sig) + 
        geom_point(shape=1) + geom_jitter() + geom_abline(colour = "grey", size = 0.5) + 
        xlim(c(-1.5, 1)) + ylim(c(-1.5, 1)) + ggtitle("Long Genes") + theme_bw() +
        geom_hline(aes(yintercept=0), linetype="dotted") + 
        geom_vline(aes(xintercept=0), linetype="dotted") +        
        geom_hline(aes(yintercept=log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=log2(1.2)), linetype="dotted") +
        geom_hline(aes(yintercept=-log2(1.2)), linetype="dotted") + 
        geom_vline(aes(xintercept=-log2(1.2)), linetype="dotted") +
        xlab(paste("LogFC RNA-Seq Values")) + ylab(paste("LogFC NanoString Values")) + 
        theme(legend.text = element_text(size = 12, face = "bold"), 
              legend.title = element_text(size = 12, colour = "black", face = "bold"),
              axis.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(size = 18, face = "bold", color = "black"),
              axis.text.y = element_text(size = 18, face = "bold", color = "black"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
              plot.title = element_blank())


long.up <- sum(res.seqnano$abs.logFC > 0 & degs.annot$long.short == "Long")
long.down <- sum(res.seqnano$abs.logFC <=  0 & degs.annot$long.short == "Long")  
short.up <- sum(res.seqnano$abs.logFC > 0 & degs.annot$long.short == "Short")
short.down <- sum(res.seqnano$abs.logFC <=  0 & degs.annot$long.short == "Short")  

ggplot(data = res.seqnano, aes(x = gene.length/1000, y = abs.logFC, col = long.short)) + 
  geom_point(size = 2) + geom_hline(aes(yintercept = 0), 
                                    colour="#FF0000", linetype="dashed", 
                                    size = 1) + 
  scale_x_continuous(trans = log10_trans(), breaks = c(1,10,100,1000)) +
  ylab(paste("abs(RNA-seq - Nanostring) LogFC (shrunken)")) + 
  xlab(paste("Gene Length")) + theme_bw() + ylim(c(-0.5,0.5)) +
  annotate("text",x = 500,y = 0.35,label = long.up, size = 7, fontface = "bold") +
  annotate("text",x = 500,y = -0.35,label = long.down, size = 7, fontface = "bold") +
  annotate("text",x = 5,y = 0.35,label = short.up, size = 7, fontface = "bold") +
  annotate("text",x = 5,y = -0.35,label = short.down, size = 7, fontface = "bold") +
  theme(legend.position="none", 
        axis.title = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 24, face = "bold", color = "black"),
        axis.text.y = element_text(size = 24, face = "bold", color = "black"))

