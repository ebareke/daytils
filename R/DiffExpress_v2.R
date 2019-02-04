#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2") ; library(DESeq2)
#biocLite("ggplot2") ; library(ggplot2)
#biocLite("clusterProfiler") ; library(clusterProfiler)
#biocLite("biomaRt") ; library(biomaRt)
#biocLite("ReactomePA") ; library(ReactomePA)
#biocLite("DOSE") ; library(DOSE)
#biocLite("KEGG.db") ; library(KEGG.db)
#biocLite("org.Mm.eg.db") ; library(org.Mm.eg.db)
#biocLite("org.Hs.eg.db") ; library(org.Hs.eg.db)
#biocLite("pheatmap") ; library(pheatmap)
#biocLite("genefilter") ; library(genefilter)
#biocLite("RColorBrewer") ; library(RColorBrewer)
#biocLite("GO.db") ; library(GO.db)
#biocLite("topGO") ; library(topGO)
#biocLite("dplyr") ; library(dplyr)
#biocLite("gage") ; library(gage)
#biocLite("ggsci") ; library(ggsci)
#########################################################################################
### If everything is installed ... ################
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
#########################################################################################

# Import gene counts table
# - skip first row (general command info) - if counts are from featureCounts
# - make row names the gene identifiers
rawdata <- read.table("deseq2_counts.txt", header = TRUE, skip = 0, row.names = 1)
countdata <- rawdata[rowSums(rawdata)>10, ]
head(countdata)

# Import metadata file
# - make row names the matching sampleID's from the countdata
metadata <- read.delim("metadata.txt", sep = ",", row.names = 1)

# Add sampleID's to the mapping file
metadata$sampleid <- row.names(metadata)

# Reorder sampleID's to match featureCounts column order. 
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]

# Make sure ID's are correct
head(metadata)

# - countData : count dataframe
# - colData : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare
ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~Group)

# Find differential expressed genes
ddsMat <- DESeq(ddsMat)

# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, contrast=c("Group", "HOM", "HET"), pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)

## 
# out of 25097 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1602, 6.4%
# LFC < 0 (down)     : 678, 2.7%
# outliers [1]       : 30, 0.12%
# low counts [2]     : 9612, 38%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Check directionality of the log2 fold changes
## Log2 fold change is set as (mutantMm / wildType)
## Positive fold changes = Increased in mutantMm
## Negative fold changes = Decreased in mutantMm
mcols(results, use.names = T)

# DataFrame with 6 rows and 2 columns
                       type                                        description
                <character>                                        <character>
# baseMean       intermediate          mean of normalized counts for all samples
# log2FoldChange      results log2 fold change (MLE): Group mutantMm vs wildType
# lfcSE               results         standard error: Group mutantMm vs wildType
# stat                results         Wald statistic: Group mutantMm vs wildType
# pvalue              results      Wald test p-value: Group mutantMm vs wildType
# padj                results                              fdr adjusted p-values


# Mouse genome database (Select the correct one)
library(org.Mm.eg.db) 

# Add gene full name
results$description <- mapIds(x = org.Mm.eg.db, keys = row.names(results), column = "GENENAME", keytype = "SYMBOL", multiVals = "first")

# Add gene symbol
results$symbol <- row.names(results)

# Add ENTREZ ID
results$entrez <- mapIds(x = org.Mm.eg.db, keys = row.names(results), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Add ENSEMBL
results$ensembl <- mapIds(x = org.Mm.eg.db, keys = row.names(results), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")

# Subset for only significant genes (p < 0.1)
results_sig <- subset(results, padj < 0.1)
head(results_sig)


# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), file = 'normalized_counts.txt', sep = '\t', quote = F, col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), file = 'normalized_counts_significant.txt', sep = '\t', quote = F, col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(results), file = "results_gene_annotated.txt", sep = '\t', quote = F, col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), file = "results_gene_annotated_significant.txt", sep = '\t', quote = F, col.names = NA)
            
# Write the counts and annotated results table to a .txt file
result_combined <- merge(as.data.frame(counts(ddsMat, normalized=FALSE)), as.data.frame(results), by="row.names", sort=FALSE)
write.table(x = as.data.frame(result_combined), file = "results_diff_gene_expression_analysis.txt", sep = '\t', quote = F, col.names = NA)



# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
png("diffexpr-PCAplot.png")
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 5) + # Increase point size
  scale_y_continuous(limits = c(-15, 15)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", subtitle = "Top 500 most variable genes")
dev.off()

# Heatmap
# Gather 50 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:50, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group), 
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(HOM = "lightblue", HET = "darkorange"),
  Replicate = c(Rep1 = "darkred", Rep2 = "forestgreen", Rep3 = "grey")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
png("diffexpr-Heatmap.png")
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 55, # Make the cells wider
         show_colnames = F)
dev.off()



# Volcanoplot

# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05) or -log10 (1 = 0.1)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "NonSignificant"))
                                       
# Write the counts and simplified results table to a .txt file
data2 = data
rownames(data2) = data2$gene
result_simplified <- merge(as.data.frame(counts(ddsMat, normalized=FALSE)), as.data.frame(data2[,c(2:4)]), by="row.names", sort=FALSE)
write.table(x = as.data.frame(result_simplified), file = "results_simplified_de_report.txt", sep = '\t', quote = F, col.names = NA)

# Make a basic ggplot2 object with x-y values
png("diffexpr-Volcanoplot.png")
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality", values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("HOM" / "HET"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  
print(vol)

dev.off()

#ggsave("diffexpr-Volcanoplot.png", plot = last_plot(), width = 20, height = 20, units = "cm")   



# MA plot
png("diffexpr-MAplot.png")
plotMA(results, ylim = c(-5, 5))
dev.off()


# Dispersion plot
png("diffexpr-Dispersionplot.png")
plotDispEsts(ddsMat)
dev.off()


### Finding pathways
# Set up matrix to take into account EntrezID's and fold changes for each gene
# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig_entrez$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)


## Enrich genes using the KEGG database

kegg_enrich <- enrichKEGG(gene = names(gene_matrix), organism = 'mouse', pvalueCutoff = 0.1, qvalueCutoff = 0.10)

# Plot results

png("KEGG_Enrichment_Pathways.png")
barplot(kegg_enrich, drop = TRUE, showCategory = 10, title = "KEGG Enrichment Pathways", font.size = 10)
dev.off()



## Enrich genes using the Gene Onotlogy

go_enrich <- enrichGO(gene = names(gene_matrix), OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = 0.1, qvalueCutoff = 0.10)

# Plot results
png("GO_Biological_Pathways.png")
barplot(go_enrich, drop = TRUE, showCategory = 10, title = "GO Biological Pathways", font.size = 10)
dev.off()
