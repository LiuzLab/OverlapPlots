## Calculating the gene length from the file: refGene.txt
## URL: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/
## Reason: ~350 genes were not included in the RefSeq.mm10.Gene.TransciptLength.txt file used
## in the previous analysis 

## libraries
rm(list = ls())
library(annotate)
library(GenomicRanges)

## datasets  -- just a check
dat.ncbi <- read.table(file = "../dat-info/ncbiRefSeq.txt", header = T, sep = "\t",
                  stringsAsFactors = F, quote = "")
dat.ucsc <- read.table(file = "../dat-info/refGene.txt", header = T, sep = "\t",
                       stringsAsFactors = F, quote = "")
dat.previous <- read.table(file = "../dat-info/RefSeq.mm10.Gene.TransciptLength_2018manuscript.txt", 
                           header = T, sep = "\t", stringsAsFactors = F, quote = "")
boxer.dat <- read.table("../dat/publication/GSE128178_10WT_10MeCP2_KO_whole_cell_RNAseq_exon_counts.txt",
                  sep = "\t", stringsAsFactors=F, header = TRUE, row.names = 1)
genes.boxer <- rownames(boxer.dat)
genes.boxer[!genes.boxer %in% unique(dat.ncbi$name2)] ## all are present 
genes.boxer[!genes.boxer %in% unique(dat.ucsc$name2)] ## 865 genes missing
genes.boxer[!genes.boxer %in% unique(dat.previous$gene.name)] ## 358 genes missing

## compute gene length
dat <- dat.ncbi
rm(dat.previous, dat.ucsc)
dim(dat)
# chrs <- unique(dat$chrom)
# chrs <- chrs[grep(pattern = "chr([0-9]+|[X|Y|M])$", x = chrs)]
# dat <- dat[dat$chrom %in% chrs,]
dim(dat)

## removing .(.) in NM ids
## now calculating the gene length per NM id
dat$name <- gsub(pattern = "(.*)\\.(.*)", replacement = "\\1", x = dat$name)
refseq <- data.frame(data.frame(matrix(ncol = 6, nrow = 0)))
colnames(refseq) <- c("chr", "tx.start", "tx.end", "strand",
                      "gene.name","gene.length")
genes <- sort(unique(dat$name2))
for(i in 1:length(genes)){
  gene <- genes[i]
  cat("processing",i,"for gene:",gene,"\n")
  chr <- unique(dat$chrom[dat$name2 %in% gene])
  str <- unique(dat$strand[dat$name2 %in% gene])
  if(length(chr) == 1 & length(str) == 1){
    refseq[i,"chr"] <- chr
    refseq[i,"tx.start"] <- min(dat[dat$name2 == gene, "txStart"])
    refseq[i,"tx.end"] <- max(dat[dat$name2 == gene, "txEnd"])
    refseq[i,"strand"] <- str
    refseq[i,"gene.name"] <- gene
    refseq[i,"gene.length"] <- refseq[i,"tx.end"]-refseq[i,"tx.start"]
  }
  if(length(chr) == 1 & length(str) > 1){
    mat <- dat[dat$name2 == gene,]
    mat1 <- mat[mat$strand == "+", ]
    mat2 <- mat[mat$strand == "-", ]
    tx.start.pos <- min(mat1[,"txStart"])
    tx.end.pos <- max(mat1[,"txEnd"])
    tx.start.neg <- min(mat2[,"txStart"])
    tx.end.neg <- max(mat2[,"txEnd"])
    gl.pos <- tx.end.pos - tx.start.pos
    gl.neg <- tx.end.neg - tx.start.neg
    if(gl.pos > gl.neg){
      refseq[i,"chr"] <- unique(mat1$chrom)
      refseq[i,"tx.start"] <- tx.start.pos
      refseq[i,"tx.end"] <- tx.end.pos
      refseq[i,"strand"] <- "+"
      refseq[i,"gene.name"] <- gene
      refseq[i,"gene.length"] <- gl.pos
    }else{
      refseq[i,"chr"] <- unique(mat1$chrom)
      refseq[i,"tx.start"] <- tx.start.neg
      refseq[i,"tx.end"] <- tx.end.neg
      refseq[i,"strand"] <- "-"
      refseq[i,"gene.name"] <- gene
      refseq[i,"gene.length"] <- gl.neg
    }
    rm(mat,mat1,mat2,tx.start.pos,tx.end.pos,tx.start.neg,
       tx.end.neg,gl.pos,gl.neg)
  }
  if(length(chr) > 1 & length(str) == 1){ ## there are 2 genes: G530011O06Rik, Erdr1 
    mat <- dat[dat$name2 == gene,]
    mat$gene.length <- mat$txEnd - mat$txStart  
    mat <- mat[which(mat$gene.length == max(mat$gene.length)),]
    refseq[i,"chr"] <- mat$chrom
    refseq[i,"tx.start"] <- mat$txStart
    refseq[i,"tx.end"] <- mat$txEnd
    refseq[i,"strand"] <- mat$str
    refseq[i,"gene.name"] <- gene
    refseq[i,"gene.length"] <- refseq[i,"tx.end"]-refseq[i,"tx.start"]
  }
  if(length(chr) > 1 & length(str) > 1){
    cat("look at this gene:",gene,"\n")
    break
  }  
  rm(gene, str)
}

## checks
sum(is.na(refseq$gene.length))
genes.boxer[!genes.boxer %in% refseq$gene.name]

## writing it
write.table(x = refseq, file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)
save(refseq, file = "../dat-info/refseq_ncbi_boxerdataset/mm10_ncbi-refSeqGene_Dec2019.RData")

## comparing with previous version
## and compared with NCBI list and it looks good!

## RefSeq.mm10.Gene.TransciptLength_2018manuscript.txt
refseq.annot <- read.table(file = "../dat-info/refseq_ucsc_dec2019/mm10_refSeqGene_As-per-Dec2019.txt", 
                           sep = "\t", header = T, stringsAsFactors = F)
refseq.annot <- refseq.annot[,c(1:6)]
refseq.annot <- unique(refseq.annot)
genes.annot <- unique(refseq.annot$gene.name[duplicated(refseq.annot$gene.name)])
idx <- c()
for(i in 1:length(genes.annot)){
  gene <- genes.annot[i]
  mat <- refseq.annot[refseq.annot$gene.name == gene,]
  rm.ind <- rownames(mat[which(mat$gene.length != max(mat$gene.length)),])
  idx <- c(idx,rm.ind)
}
refseq.annot <- refseq.annot[which(!rownames(refseq.annot) %in% idx),]
refseq1 <- inner_join(refseq, unique(refseq.annot), by="gene.name")
refseq1$gene.length.diff <- refseq1$gene.length.x-refseq1$gene.length.y

