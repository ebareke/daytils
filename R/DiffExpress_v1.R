#
R
countData <- read.table("deseq2.counts", header=TRUE, row.names=1)
#quant <- apply(countdata,1,quantile,0.75)
#keep <- which((quant >= 10) == 1)
#countdata <- countdata[keep,]
countData <- countData[rowSums(countData)>10, ]
colData = read.csv("metadata.csv", sep="\t", row.names=1)
###
library(DESeq2)
##
dds = DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds = DESeq(dds)
#
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()
#
rld <- rlogTransformation(dds)
library( "genefilter" )
library( "gplots" )
library( "RColorBrewer" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
png("qc-heatmap2-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[colData(rld)$condition ] )
dev.off()
##
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(colData$condition))])
##
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none", col=colorpanel(100, "black", "white"),ColSideColors=mycols[colData$condition], RowSideColors=mycols[colData$condition], margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
##
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
#
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()
##
# Get differential expression results
res <- results(dds,contrast=c("condition", "mutantMm", "wildType"))
table(res$padj<0.1)
## Order by adjusted p-value
res <- res[order(res$padj), ]
library("AnnotationDbi")
library("org.mm.eg.db")
##res$symbol <- mapIds(org.Mm.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
##res$entrez <- mapIds(org.Mm.eg.db, keys=row.names(res), column="ENTREZID",keytype="ENSEMBL", multiVals="first")
##res$entrez <- mapIds(EnsDb.Mmusculus.v75, keys=row.names(res), column="ENTREZID",keytype="ENSEMBL", multiVals="first")
resOrdered <- res[order(res$padj),]
summary(resOrdered)
## Merge with normalized count data
resdata <- merge(as.data.frame(counts(dds, normalized=FALSE)), as.data.frame(res), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resOrderedData <- merge(as.data.frame(counts(dds, normalized=FALSE)), as.data.frame(resOrdered), by="row.names", sort=FALSE)
names(resOrderedData)[1] <- "Gene"
#head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")
write.csv(resOrderedData, file="diffexpr-resultsOrdered.csv")

##
maplot <- function (res, thresh=0.1, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1000, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot2.png", 1000, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=2, sigthresh=0.1, textcx=.8, xlim=c(-10, 10))
dev.off()
