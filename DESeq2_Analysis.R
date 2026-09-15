# Loading packages
library(tidyverse)
library(stringr)
library(dplyr)
library(DESeq2)           
library(ashr)
library(pheatmap)          
library(RColorBrewer)      
library(ggplot2)           
library(EnhancedVolcano)   
library(clusterProfiler)   
library(org.Hs.eg.db)      
library(AnnotationDbi)     

# ---------------------------------------------------------------------

# Reading featureCounts files
# construct vector of all file paths
count_files <- list.files(
  path = "H:/My Drive/jroverflow/retOrgs/featureCounts",
  pattern = "*.tabular.txt",
  full.names = TRUE
)
count_files <- str_sort(count_files, numeric = TRUE)
# similar vector of just sample names
count_names <- list.files(
  path = "H:/My Drive/jroverflow/retOrgs/featureCounts",
  pattern = "*.tabular.txt",
  full.names = FALSE
)
count_names <- str_split_i(count_names, ".tabular.txt", 1)
count_names <- str_sort(count_names, numeric = TRUE)

# first file as leftmost col
count_matrix <- read.table(count_files[1], TRUE, "\t")
# extract name from first file
first_sample <- count_names[1]
# renaming counts col
value_col <- setdiff(names(count_matrix), "Geneid")
names(count_matrix)[names(count_matrix) == value_col] <- first_sample

# loop, since first file already done, start at file 2
for (i in 2:length(count_files)) {
  # read file
  sample <- read.table(count_files[i], TRUE, "\t")
  # get sample name
  name <- count_names[i]
  value_col <- setdiff(names(sample), "Geneid")
  names(sample)[names(sample) == value_col] <- name
  # left join by Geneid
  count_matrix <- left_join(count_matrix, sample, "Geneid")
}
#-----------------------------------------------------------------------

# Step 4: Filtering
# Condition: Counts > 0
filtered_count_matrix <- count_matrix[rowSums(count_matrix[, -1]) > 0, ]
# omitting Geneid col as a valid col
rownames(filtered_count_matrix) <- filtered_count_matrix$Geneid
filtered_count_matrix$Geneid <- NULL

#-----------------------------------------------------------------------
# Step 5: Building DESeq2 Object
## 3 required items: count matrix, metadata, and design condition
# manually setting up metadata df
sample_timepoints <- c(
                  "D00", "D00", "D00",
                  "D10", "D10", "D10", "D10",
                  "D25", "D25", "D25", "D25",
                  "D65", "D65", "D65", "D65",
                  "D100", "D100", "D100", "D100",
                  "D180", "D180", "D180", "D180",
                  "D180", "D180",
                  "D280", "D280", "D280")
metadata <- data.frame(Sample = count_names, Timepoint = sample_timepoints)
rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL

# must be true! uncomment to check
# print(colnames(filtered_count_matrix) == rownames(metadata))

# building object
dds <- DESeqDataSetFromMatrix(filtered_count_matrix,
                              colData = metadata,
                              design = ~ Timepoint)


#------------------------------------------------------------------------
### To replicate the paper's findings, we'll contrast D25 v D65
dds <- DESeq(dds)
dds <- dds[rowSums(counts(dds)) > 0, ]
results <- results(dds, contrast = c("Timepoint", "D65", "D25"))
results <- results[order(results$padj), ]

#-------------------------------------------------------------------------
# Step 7: Analysis
# adding significance column
results$significant <- results$padj < 0.05 & abs(results$log2FoldChange) > 0.5

# annotating GeneIds
ids <- mapIds(org.Hs.eg.db, 
              keys = rownames(results),
              column = "SYMBOL",
              keytype = "ENSEMBL",
              multiVals = "first")
results$symbol <- ids
results$symbol[is.na(results$symbol)] <-
  rownames(results)[is.na(results$symbol)]
rownames(results) <- results$symbol
results$symbol <- NULL
head(results)


# normalized counts
# note: this is a matrix, not a df,  all values in the matrix are numerical
normalized_counts <- counts(dds, normalized = TRUE)
rownames(normalized_counts) <- mapIds(org.Hs.eg.db,
                                      keys = rownames(normalized_counts),
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")

# ----------------------------------------------------------------------
# PCA plot
vsd <- vst(dds)

png("H:/My Drive/jroverflow/retOrgs/toGitHub/PCAPlot.png", width = 2000,
    height = 2000,
    res = 300)

pca_plot <- plotPCA(vsd, intgroup = "Timepoint")
pca_plot

dev.off()

#------------------------------------------------------------------------
# Sample-to-sample distance heatmap
rownames(vsd) <- ids
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)

png("H:/My Drive/jroverflow/retOrgs/toGitHub/sts-heatmap.png", 
    width =  3000,    
    height = 3000,
    res = 300)
sts_heatmap <- pheatmap(sampleDistMatrix,
                        clustering_distance_rows = sampleDists,
                        clustering_distance_cols = sampleDists,
                        labels_row = rownames(sampleDists),
                        labels_col = colnames(sampleDists),
                        col = colors)
dev.off()

#-------------------------------------------------------------------------
# Heat map: version 1, sample by sample
# selects top 50 (numerically) expressed genes
select <- order(rowMeans(normalized_counts), decreasing = TRUE)[1:50]

# essentially the metadata (reproduced for pheatmap)
df <- as.data.frame(colData(dds)[, c("Timepoint")])
df <- data.frame(Timepoint = dds$Timepoint)
rownames(df) <- colnames(normalized_counts)


heat_map <- pheatmap(normalized_counts[select,], cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE, annotation_col = df)
heat_map

#-----------------------------------------------------------------------------
# volcano plot
res_shrunk <- lfcShrink(dds, 
                        contrast = c("Timepoint", "D25", "D65"), 
                        type = "ashr")

# annotating GeneIds
res_shrunk$symbol <- mapIds(org.Hs.eg.db, 
                            keys = rownames(res_shrunk),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
res_shrunk$symbol[is.na(res_shrunk$symbol)] <- rownames(res_shrunk)[is.na(res_shrunk$symbol)]
rownames(res_shrunk) <- res_shrunk$symbol
res_shrunk$symbol <- NULL
# similar labels used by authors
to_label <- c("LHX4","ARR3","THRB","CRX","VSX2","POU4F2","RXRG","ISL1","ATOH7",
                              "HES1","GNGT2","SIX6","GMNN","NR2E1","DCX","NES","MKI67","PAX6",
                              "SIX3","GNL3")
# also labelling top 20 most sig DE genes to label
top <- rownames(res_shrunk[order(res_shrunk$padj, decreasing = FALSE),])[1:20]
to_label <- c(to_label, top)

# colors
colors <- c("grey", "#bdc9e1","#74a9cf", "#103F5C")

png("H:/My Drive/jroverflow/retOrgs/toGitHub/volcanoPlot.png", width = 3000,
    height = 3000,
    res = 300)
EnhancedVolcano(res_shrunk,
                lab = rownames(res_shrunk),
                x = 'log2FoldChange',
                y = 'padj',
                title = "D65 vs D25",
                subtitle = NULL,
                selectLab = to_label,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 3.0,
                col = colors,
                colAlpha = 0.85,
                boxedLabels = TRUE,
                drawConnectors = TRUE)


dev.off()


#-------------------------------------------------------------------------
# Gene Ontology (GO) enrichment analysis
# ordering genes by lowest padj value (eg most statistically sig)
res_ordered <- res_shrunk[order(res_shrunk$padj),]
res_df <- as.data.frame(res_ordered)
res_df$ensembl <- rownames(res_df)
# annotating entrezid
res_df$entrez <- mapIds(org.Hs.eg.db, 
              keys = rownames(res_df),
              column = "ENTREZID",
              keytype = "SYMBOL",
              multiVals = "first")

sig_entrez <- na.omit(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 0.5)$entrez)
# na.omit removes "NA" values present in the matrix
universe_entrez <- na.omit(res_df$entrez)

# running GO enrichment: ont can be either
# 1. BP = biological processes (most common)
# 2. MF = molecular function
# 3. CC = cellular component
png("H:/My Drive/jroverflow/retOrgs/toGitHub/goEnrichmentPlot.png", 
    width = 3000,
    height = 3000,
    res = 300)
# gene parameter: significantly expressed DE genes
# universe parameter: umbrella of all genes tested for DE (eg count matrix)
# pAdjustMethod will almost always be BH (Benjamini-Hochberg)
go_results <- enrichGO(gene = sig_entrez,
                        universe = universe_entrez,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        readable = TRUE)

go_plot <- dotplot(go_results, showCategory = 15)
go_plot
 
dev.off()

#------------------------------------------------------------------------
# Kegg enrichment analysis

png("H:/My Drive/jroverflow/retOrgs/toGitHub/keggEnrichmentPlot.png", 
    width = 3000,
    height = 3000,
    res = 300)
kegg_results <- enrichKEGG(gene = sig_entrez,
                            universe = universe_entrez,
                            organism = "hsa",
                            pvalueCutoff = 0.05)

kegg_plot <- dotplot(kegg_results, showCategory = 15)
kegg_plot

dev.off()

#----------------------------------------------------------------------------
# GSEA 
# since mapIds() didn't label everything, reextracting results
re_results <- results(dds, contrast = c("Timepoint", "D65", "D25")) 
# ordering is done either by log2fc or stat, stat little more robust
res_ranked <- re_results[order(re_results$stat, decreasing = TRUE),]
# pairing stat nums w associated genes
ranked_genes <- res_ranked$stat
names(ranked_genes) <- rownames(res_ranked)
# running gseGO
gsea_results <- gseGO(ranked_genes,
                      ont = "BP",
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      pvalueCutoff = 0.05,
                      eps = 0,
                      verbose = TRUE)
as.data.frame(gsea_results) 

# select top 10 upreg pathways by nes
top_upreg <- gsea_results[order(gsea_results$NES, decreasing = TRUE)[1:10], ]
# same with downreg pathways
top_downreg <- gsea_results[order(-gsea_results$NES, decreasing = TRUE)[1:10],]
# combine two to one df
top_pathways <- rbind(top_upreg, top_downreg) |> select_()
as.data.frame(top_pathways)
# simplify table for graphability
top_pathways |> dplyr::select(ID, Description, NES, p.adjust, qvalue)
# ensure ggplot won't alphabetically sort descriptions
top_pathways$Description <- factor(top_pathways$Description, levels = top_pathways$Description)

png("H:/My Drive/jroverflow/retOrgs/toGitHub/gsea_top20pathways.png", 
    width = 3000,
    height = 1500,
    res = 300)

ggplot(top_pathways, aes(x = NES, y = Description, color = NES)) +
  geom_segment(aes(x = 0, xend = NES, y = Description, linewidth = 5)) + 
  scale_color_gradient2(low = "#af8dc3",
                        mid = "#f7f7f7",
                        high = "#7fbf7b") + 
  labs(x = "NES", y = "Pathway Description") + theme_minimal()

dev.off()
