############ Make a Heatmap (based on the tutorial here https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/)
# Load libraries
library(readxl);library(DESeq2);library(reshape2);library(ggdendro);library(khroma);library(viridis);library(gridExtra);library(cowplot);library(pals)

# Perform variance stabilizing transformation
vsd <- vst(ddsObj)
vst_mat <- assay(vsd)

# Get tissue assignment for each sample
tissue_info <- colData(ddsObj)$Tissue
names(tissue_info) <- colnames(vst_mat)  # label with sample names

# Create a list mapping each tissue to the sample columns
tissue_groups <- split(names(tissue_info), tissue_info)

# Collapse vst matrix by computing rowMeans per tissue
collapsed_vst <- sapply(tissue_groups, function(samples) {
  rowMeans(vst_mat[, samples, drop = FALSE], na.rm = TRUE)
})

# Convert to data frame
collapsed_vst_df <- as.data.frame(collapsed_vst) %>%
  rownames_to_column("Geneid")


# Set rownames and remove Geneid column before scaling
mat <- collapsed_vst_df
rownames(mat) <- mat$Geneid
mat$Geneid <- NULL

# Z-score across tissues (i.e., scale by row)
Z <- t(scale(t(as.matrix(mat))))


# Convert to dataframe and keep gene IDs
Z <- as.data.frame(Z) %>%
  rownames_to_column("Geneid")

normcounts <- read_excel("Plodia_SG_SalG_Diff_Exp_normalized_counts.xlsx", sheet = "MSG_vs_PSG_diff_exp_4tissue") 

normcounts <- normcounts %>%
  select(Geneid, `max(msg,psg)`)

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(normcounts, by = c("Geneid" = "Geneid"))


# Filter significant genes (adjust as needed)
sigGenes <- data.frame(Geneid = topGenesMSG_PSG$Geneid[
  topGenesMSG_PSG$padj <= 0.01 & abs(topGenesMSG_PSG$log2FoldChange) > 3 & topGenesMSG_PSG$`max(msg,psg)` > 300
])



Z <- Z[Z$Geneid %in% sigGenes$Geneid, ]

flybasenames <- topGenesMSG_PSG %>%
  select(Geneid, FB_gene_symbol) %>%
  drop_na() %>%
  add_count(Geneid) %>%
  filter(n == 1) %>%
  select(-n)

flybasenames_unique <- flybasenames %>%
  group_by(FB_gene_symbol) %>%
  filter(n() == 1) %>%     # Keep only those with a single occurrence
  ungroup()



Z <- Z %>%
  left_join(flybasenames_unique, by = c("Geneid" = "Geneid")) %>%
  mutate(
    Geneid = ifelse(!is.na(FB_gene_symbol), FB_gene_symbol, Geneid)  # Replace Geneid with FB_symbol if match exists
  ) %>%
  select(-FB_gene_symbol) 



# Melt to long format
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")


# Convert to wide matrix format for clustering
Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL


# Compute distances and clusters
distanceGene <- dist(Z_df_matrix)
distanceSample <- dist(t(Z_df_matrix))
clusterSample <- hclust(distanceSample, method = "average")
clusterGene <- hclust(distanceGene, method = "average")

# Construct a dendogram for samples
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

#make dendogram for genes 
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
geneDendrogram <- ggplot(geneDendrogramData) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dendro()

# Re-factor samples for ggplot2 for ordering 
Z_df$Sample  <- factor(Z_df$Sample , levels=c("MSG","PSG","SalG","Head"))
Z_df$Gene <- factor(Z_df$Gene, levels=clusterGene$labels[clusterGene$order])

#define your color palette 
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm

# Now add your manual symbol annotations, download googlesheet as excel sheet

manual_annotations <- read_excel("ilPloInte3.2_manually_curated_gene_table.xlsx")

manual_annotations <- manual_annotations %>%
  select(`NCBI Gene ID`, Symbol)

# Preserve factor levels from clustering
gene_levels <- levels(Z_df$Gene)

# Create mapping vector for label replacement
label_map <- manual_annotations %>%
  filter(`NCBI Gene ID` %in% gene_levels) %>%
  distinct(`NCBI Gene ID`, Symbol) %>%
  deframe()

# Apply label map to factor levels (not values)
new_levels <- ifelse(gene_levels %in% names(label_map), label_map[gene_levels], gene_levels)
levels(Z_df$Gene) <- new_levels

# Construct the heatmap
heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))

# Clean up the plot layout with cowplot 
# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = coolwarm(256)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"
    ))

heatmap
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots
#top_row <- plot_grid(NULL, sampleDendrogram, ncol = 2, rel_widths = c(0.2, 1))
mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))

write.table(Z_df_matrix, file="heatmap_raw_table_230genes.tsv", quote=F, sep="\t",row.names=TRUE, na="")
