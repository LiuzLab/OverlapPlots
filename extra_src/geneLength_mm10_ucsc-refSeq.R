## Calculating the gene length from the file: refGene.txt
## URL: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/
## Reason: ~350 genes were not included in the RefSeq.mm10.Gene.TransciptLength.txt 
## file used in previous analysis

## libraries
rm(list = ls())
library(annotate)
library(GenomicRanges)

## dataset
dat <- read.table(file = "../dat-info/refGene.txt", header = T, sep = "\t", 
                  stringsAsFactors = F, quote = "")
dim(dat)
chrs <- sort(unique(dat$chrom)[1:21])
dat <- dat[dat$chrom %in% chrs,]
dim(dat)

## some checks; check name/Refseq symbol
refseq.nm <- unique(dat[duplicated(dat$name), "name"])
dat1 <- dat[dat$name %in% refseq.nm, ]
dim(dat1)
idx.rm <- c()
for(i in 1:length(refseq.nm)){
  mat <- dat1[dat1$name %in% refseq.nm[i],] 
  mat$gene.length <- mat$txEnd-mat$txStart
  chr.same <- length(unique(mat$chrom)) == 1
  strand.same <- length(unique(mat$strand)) == 1
  gl.same <- length(unique(mat$gene.length)) == 1
  if(gl.same == TRUE){
    idx <- rownames(mat)[2:nrow(mat)]
    idx.rm <- c(idx.rm, idx) 
  }
  if(gl.same == FALSE){
    ## keeping the biggest gene length for that NM (and the first one only)
    idx1 <- which(mat$gene.length != max(mat$gene.length))
    idx2 <- which(mat$gene.length == max(mat$gene.length))                   
    if(length(idx2) > 1){
      idx <- c(idx1,idx2[2:length(idx2)])
      idx <- rownames(mat[idx,])  
    }
    else{
      idx <- rownames(mat[which(mat$gene.length != max(mat$gene.length)),])
    }
    idx.rm <- c(idx.rm, idx)
  }
}
rm(mat, idx, chr.same, strand.same, gl.same, i)
dat <- dat[which(!rownames(dat) %in% idx.rm),]

## now calculating the gene length per NM id
refseq <- data.frame(data.frame(matrix(ncol = 6, nrow = 0)))
colnames(refseq) <- c("chr", "tx.start", "tx.end", "strand",
                      "gene.name","gene.length")
genes <- unique(dat$name2)
for(i in 1:length(genes)){
  gene <- genes[i] ## genes[genes %in% "Gm3752"]
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
  else if(length(str) > 1){
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
      refseq[i,"gene.name"] <- unique(mat1$name2)
      refseq[i,"gene.length"] <- gl.pos
    }else{
      refseq[i,"chr"] <- unique(mat1$chrom)
      refseq[i,"tx.start"] <- tx.start.neg
      refseq[i,"tx.end"] <- tx.end.neg
      refseq[i,"strand"] <- "-"
      refseq[i,"gene.name"] <- unique(mat1$name2)
      refseq[i,"gene.length"] <- gl.neg
    }
  }
  
  # refseq[i,"chr"] <- dat[i, "chrom"]
  # refseq[i,"tx.start"] <- max(dat[dat$name2 == gene, "txStart"])
  # refseq[i,"tx.end"] <- max(dat[dat$name2 == gene, "txEnd"])
  # refseq[i,"strand"] <- unique(dat[dat$name2 == gene, "strand"])
  # refseq[i,"gene.name"] <- gene
  # refseq[i,"gene.length"] <- refseq[i,"tx.end"]-refseq[i,"tx.start"]
  # 
  # ## check if chr and strand is unique
  # if(length(refseq$chr[i]) > 1){
  #   cat("Error: multiple chr for a gene \n")
  #   break;
  # }
  
  # ## now calculating the max transcript.length (have to loop per transcript and not per gene)
  # transcripts <- unique(dat[dat$name2 == gene, "name"])
  # exon.starts <- strsplit(dat[dat$name2 == gene, "exonStarts"], ",")
  # exon.ends <- strsplit(dat[dat$name2 == gene, "exonEnds"], ",")
  # for(j in 1:length(transcripts)){
  #   refseq[i,"transcript"] <- transcripts[j]
  #   refseq[i,"transcript_length"] <- as.numeric(exon.ends[[j]])-as.numeric(exon.starts[[j]])
  # }
}
refseq <- refseq[!is.na(refseq$gene.length), ]
write.table(x = refseq, file = "../dat-info/refseq_ucsc_dec2019/mm10_refSeqGene_As-per-Dec2019.txt", 
            sep = "\t", quote = F, col.names = T, row.names = F)
save(refseq, file = "../dat-info/refseq_ucsc_dec2019/mm10.refSeqGene.RData")

## comparing with previous version
## and compared with UCSC list and it looks good!

refseq.annot <- read.table(file = "../dat-info/RefSeq.mm10.Gene.TransciptLength_2018manuscript.txt", 
                           sep = "\t", header = T, stringsAsFactors = F)
refseq.annot <- refseq.annot[refseq.annot$Chr_Name %in% chrs,]
refseq.annot <- refseq.annot[,c(1:6)]
refseq.annot <- unique(refseq.annot)
genes <- unique(refseq.annot$gene.name[duplicated(refseq.annot$gene.name)])
idx <- c()
for(i in 1:length(genes)){
  gene <- genes[i]
  mat <- refseq.annot[refseq.annot$gene.name == gene,]
  rm.ind <- rownames(mat[which(mat$gene.length != max(mat$gene.length)),])
  idx <- c(idx,rm.ind)
}
refseq.annot <- refseq.annot[which(!rownames(refseq.annot) %in% idx),]
refseq1 <- inner_join(refseq, unique(refseq.annot), by="gene.name")
refseq1$gene.length.diff <- refseq1$gene.length.x-refseq1$gene.length.y
