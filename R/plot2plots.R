## box plot to check if the data was normalized or not
#' Normalization check box plot
#'
#' @param data data
#' @param samples samples
#'
#' @return boxplot
#' @noRd
#'
#' @examples
boxPlot <- function(data, samples){
  par(mar=c(2+round(max(nchar(colnames(data)))/2),4,2,1),font=2)
  boxplot(data, boxwex=0.6, notch=TRUE, outline=FALSE, las=2, col=samples)
  legend("topleft", levels(samples), fill = palette(), bty="n", cex = 0.8,
         xpd = TRUE)
}

## Gene Expression Level in Microarray data
#' Gene Expression Level in Microarray data
#'
#' @param data data
#' @param gene.id gene ids
#' @param genotypes genotypes
#'
#' @return ggplot object with gene expression levels
#' @noRd
#'
#' @examples
## Gene Expression Level in Microarray data
GeneExpLevels <- function(data, gene.id, genotypes){

  ## MeCP2 expression level
  ind <- which(rownames(data) == gene.id[1])
  matGene <- as.vector(as.matrix(data[ind,]))
  plotDat <- data.frame(sampleName = colnames(data),
                        normCounts = matGene,
                        genotype = genotypes,
                        Genotype_Condition = samples)
  print(ggplot(plotDat, aes(x = sampleName, y = normCounts, fill = genotype)) +
          geom_bar(stat="identity") + ylab("Log2 Normalized Data") +
          ggtitle(paste("Barplot for Gene ",gene.id[2],sep="")) + theme_bw() +
          facet_grid(. ~ Genotype_Condition,  space = "free", scale="free") +
          theme(plot.title = element_text(size = 24, face = "bold"),
                axis.title.y= element_text(size = 20, colour = "black",
                                           face = "bold"),
                axis.text.y = element_text(size = 20, colour = "black",
                                           face = "bold"),
                axis.text.x = element_blank(), legend.position="none",
                axis.ticks.x = element_blank(), axis.title.x =  element_blank(),
                strip.text = element_text(size = 18, colour = "red",
                                          face = "bold")))
}

## MDS Plot
#' MDS Plot
#'
#' @param data data
#' @param genotypes genotypes
#' @param conditions conditions
#'
#' @return ggplot obuect MDS plot
#' @noRd
#'
#' @examples
## MDS Plot
MDSplot <- function(data, genotypes, conditions){
  mdsDist = cmdscale(d = dist(t(data)), eig = TRUE, k = 2)
  mdsDist = data.frame(genotypes, x = mdsDist$points[,1]/1e4,
                       y = mdsDist$points[,2]/1e4)

  ## plot
  if(missing(conditions)){
    print(ggplot(mdsDist, aes(x = x, y = y, color = genotypes)) +
          geom_point(size = 1) + ## shape = genotypes,
          ylab("MDS Coordinate 2 (x 1e4)") + xlab("MDS Coordinate 1 (x 1e4)") +
          theme_bw() + theme(legend.text = element_text(size = 18,
                      face = "bold"),
                      legend.title = element_text(size = 18, colour = "black",
                      face = "bold"),
                      axis.title = element_text(size = 18, face = "bold"),
                      axis.text.x = element_text(size = 18, face = "bold",
                                                          color = "black"),
                      axis.text.y = element_text(size = 18, face = "bold",
                                                          color = "black"),
                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
  else{
    print(ggplot(mdsDist, aes(x = x, y = y, shape = conditions,
                              color = genotypes)) + geom_point(size = 8) +
          ylab("MDS Coordinate 2 (x 1e4)") + xlab("MDS Coordinate 1 (x 1e4)") +
          theme_bw() + theme(legend.text = element_text(size = 18,
                                                          face = "bold"),
                      legend.title = element_text(size = 18, colour = "black",
                                                           face = "bold"),
                      axis.title = element_text(size = 18, face = "bold"),
                      axis.text.x = element_text(size = 18, face = "bold",
                                                          color = "black"),
                      axis.text.y = element_text(size = 18, face = "bold",
                                                          color = "black"),
                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
}

## PCA plot
#' PCA plot
#'
#' @param data data
#' @param genotypes genotypes
#' @param conditions conditions
#'
#' @return ggplot object PCA plot
#' @noRd
#'
#' @examples
## PCA plot
PCAplot <- function(data, genotypes, conditions){
  ## Calculating PC components
  pcs = prcomp(t(data), center = TRUE)
  percentVar = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2)
  if(missing(conditions)){
    print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = genotypes)) +
            xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) +
            geom_point(size = 8) + theme_bw() +
            theme(legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, colour = "black",
                                              face = "bold"),
                  plot.title = element_blank(),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold",
                                             color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold",
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
  else{
    print(ggplot(as.data.frame(pcs$x), aes(PC1,PC2, color = genotypes,
                                           shape = conditions)) +
            xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) +
            geom_point(size = 8) + theme_bw() +
            theme(legend.text = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size = 16, colour = "black",
                                              face = "bold"),
                  plot.title = element_blank(),
                  axis.title = element_text(size = 18, face = "bold"),
                  axis.text.x = element_text(size = 16, face = "bold",
                                             color = "black"),
                  axis.text.y = element_text(size = 16, face = "bold",
                                             color = "black"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
  }
}

#' PC variance labels
#'
#' @param x x
#' @param pc pc
#'
#' @return PC with its variance
#' @noRd
#'
#' @examples
makeLab = function(x,pc) {
  paste0("PC",pc,": ",x,"% variance")
}

## Scatter Plot for DEGs
#' Scatter Plot for DEGs
#'
#' @param dat data
#' @param log2FC log2 fold change
#' @param comp.between comp between
#' @param pval p-values
#'
#' @return scatter plot
#' @noRd
#'
#' @examples
#' \dontrun{plotScatter(dat = degs.dat[,c("gene.name","logFC", "FDR",
#' "gene.length")],
#' log2FC = log2(1), comp.between = "")}
## Scatter Plot for DEGs
plotScatter <-function(dat, log2FC, comp.between, pval = 0.05){
  colnames(dat) <- c("gene.name", "logFC", "adj.P.Val", "gene.length")
  gene.type <- ifelse((dat$adj.P.Val < pval & abs(dat$logFC) > log2FC),
                      ifelse(dat$gene.length > 100e3, "Long Genes",
                             "Short Genes"),"Not Stat. Signif.")
  ind <- which(gene.type != "Not Stat. Signif.")
  dat <- dat[ind,]
  gene.type <- gene.type[ind]

  ## Contingency Table
  up.LongGenes <- sum(dat$logFC > log2FC & dat$gene.length > 100e3)
  down.LongGenes <- sum(dat$logFC < -log2FC & dat$gene.length > 100e3)
  up.ShortGenes <- sum(dat$logFC > log2FC & dat$gene.length <= 100e3)
  down.ShortGenes <- sum(dat$logFC < -log2FC & dat$gene.length <= 100e3)
  cont.tab <- matrix(data = c(up.LongGenes, down.LongGenes, up.ShortGenes,
                              down.ShortGenes), nrow = 2, ncol = 2)
  rownames(cont.tab) <- c("Up", "Down")
  colnames(cont.tab) <- c("Long.Gene", "Short.Gene")
  print(cont.tab)

  ## qplot
  print(qplot(y = dat$logFC, x = dat$gene.length/1000, colour = gene.type,
      xlab = "Gene Length in KB", ylab = paste("Log2 Fold Change",
                                               comp.between)) +
      scale_x_continuous(trans = log10_trans(), breaks = c(0,1,10,100,1000)) +
      coord_cartesian(ylim = c(-1.5,1.5)) + theme_bw() +
      annotate("text", x = 500, y=1.5, label= cont.tab[1,1], size=7,
               fontface="bold") +
      annotate("text", x = 500, y=-1.5, label= cont.tab[2,1], size=7,
               fontface="bold") +
      annotate("text", x = 1, y=1.5, label = cont.tab[1,2], size=7,
               fontface="bold") +
      annotate("text", x = 1, y=-1.5, label = cont.tab[2,2], size=7,
               fontface="bold") +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
            legend.position="none",
            axis.text.x = element_text(size = 18, face = "bold",
                                       color = "black"),
            axis.text.y = element_text(size = 18, face = "bold",
                                       color = "black"),
            legend.text = element_text(size = 18, face = "bold"),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
}

## Scatter Plot with lm line
#' Scatter Plot with lm line
#'
#' @param dat data
#'
#' @return Scatter plot
#' @noRd
#'
#' @examples
## Scatter Plot with lm line
scatterlm <- function(dat){
  r.sq <- paste("R^2 = ",format(summary(lm(gene.length ~ logFC.crude,
                                           dat))$r.squared, digits = 2))
  print(summary(lm(gene.length ~ logFC.crude, dat)))
  cat(summary(lm(gene.length ~ logFC.crude, dat))$r.squared,"\n")
  print(qplot(y = logFC.crude, x = gene.length/1000, data = dat) +
          geom_smooth(method = "lm") + xlab("Gene Length") +
          ylab("Mean Log2 Fold Change") + theme_bw() +
          scale_x_continuous(trans = log10_trans(),
                             breaks = c(0,1,10,100,1000)) +
          annotate("text", x = 500, y=1.5, label = r.sq, size = 5,
                   fontface="bold") +
          theme(plot.title = element_text(size = 18, face = "bold"),
                axis.title.y = element_text(size = 18, colour = "black",
                                            face = "bold"),
                axis.title.x = element_text(size = 18, colour = "black",
                                            face = "bold"),
                axis.text.y = element_text(size = 18, colour = "black",
                                           face = "bold"),
                axis.text.x = element_text(size = 18, colour = "black",
                                           face = "bold"),
                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")))
}
