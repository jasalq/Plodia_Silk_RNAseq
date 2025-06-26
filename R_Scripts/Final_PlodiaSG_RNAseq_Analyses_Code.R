# Plodia Silk Gland Diff Exp Analysis (MSG vs PSG)
# Written by Jasmine D. Alqassar 2025 to use DESeq2 to analyze Plodia partitioned silk gland, salivary gland, and head RNAseq

# Set working directory
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/@LabData/2025_PlodiaSG/RNAseq_Analyses")

# Install necessary packages
install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")

# Load necessary packages
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2); library(readxl)

# Read the sample information into a data frame
sampleinfo <- read_tsv("all_samples2.txt")
sampleinfo

# Read the counts file in
seqdata <- read_tsv("Pi_SG_SalG.new.annotation.featurecounts.txt", comment="#")

# Remove columns other than GeneID and sample counts
seqdata = seqdata %>%
  select(-Chr, -Start, -End, -Strand, -Length)

# Add your curated gene symbols instead 
#first add Geneid for each old symbol and then add the new symbol
Symbol_to_Geneid <- read_excel("Plodia_gene_IDs.xlsx", col_types = "text") %>%
  select(Geneid, Symbol) %>%
  distinct(Symbol, .keep_all = TRUE)


seqdata <- seqdata %>%
  rename("Geneid" = "old_symbol") 

seqdata <- seqdata %>%
mutate(old_symbol = as.character(old_symbol))

library(readxl)
  
manual_annotations <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx")

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

seqdata <- seqdata %>%
  left_join(Symbol_to_Geneid, by = c("old_symbol" = "Symbol"))  %>%
  drop_na()


seqdata <- seqdata %>%
left_join(manual_annotations, by = c("Geneid" = "Geneid"))  %>%
  select(-`old_symbol`) 

  
  
# Simplify to sample name for each counts column 
seqdata <- seqdata %>%
  rename_with(~ ifelse(
    . %in% c("Geneid", "Symbol"),
    .,
    str_extract(., "[^/]+(?=_pass2_mappedAligned\\.sortedByCoord\\.out\\.bam)")
  ))
# Transform raw data to matrix of counts:
countdata <- seqdata  %>%
  group_by(Geneid) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# Remove the genes that are not/lowly expressed:
keep <- rowSums(countdata) > 0 
head(keep)
table(keep)

countdata <-countdata[keep, ]
dim(countdata)
head(countdata)

## QC of counts; count distributions and visualize
summary(countdata)
boxplot(countdata, las=2)

# Look at variance-mean relationship
plot(x=rowMeans(countdata), y=rowSds(countdata),
     main="sd vs mean",
     xlim=c(0,10000),
     ylim=c(0, 5000))
abline(lm(rowSds(countdata) ~ rowMeans(countdata)), col ="magenta")

############# DESeq2 Analysis #####################
####### DESeq2 Analysis MSG vs PSG

`%out%` <- function(x, y) !(x %in% y)

# remove the concatenated read samples from samples and count data
sampleinfo <- sampleinfo %>%
  filter(.[[1]] %out% c("1B_cat", "2B_cat", "3B_cat", "4B_cat"))

countdata <- countdata[, !colnames(countdata) %in% c("1B_cat", "2B_cat", "3B_cat", "4B_cat")]

# adding identifier to collapse by sequencing runs 

for(i in 1:nrow(sampleinfo)){
  sampleinfo$seq_run[i] <- strsplit(sampleinfo$Sample[i], '_')[[1]][1]}

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sampleinfo,
  design = ~ Tissue
)

ddscoll <- collapseReplicates(dds, dds$seq_run, renameCols=TRUE) #collapsing sequencing runs of PSG

ddsObj.raw <- DESeq(ddscoll)
ddsObj <- estimateSizeFactors(ddsObj.raw)

colData(ddsObj.raw)
colData(ddsObj)

ddsObj <- DESeq(ddsObj.raw)

res <- results(ddsObj, alpha=0.05) #adding p-value cutoff
res

# Display contrast names
resultsNames(ddsObj)

#Get right contrasts, in this case MSG vs PSG
Tissue_MSG_vs_PSG <- results(ddsObj, alpha = 0.05, contrast=c("Tissue","MSG","PSG"))

#Display most diff expressed genes for this contrast

sum(Tissue_MSG_vs_PSG$padj < 0.05, na.rm = TRUE)

topGenesMSG_PSG <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>% 
  head(2914)


allGenesMSG_PSG_dds <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) 

#A quick Volcano Plot 
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))

threshold <- 1000 #assigning an arbitrary cutoff to label all points above this -log10(padj) value because the plot gets too crowded

vol_plot <- topGenesMSG_PSG %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = direction)) +  
  geom_point() +  
  scale_color_manual(values = c("down" = "#26b3ff", "ns" = "grey", "up" = "#bb0c00"),  
                     labels = c("down" = "PSG", "ns" = "Not significant", "up" = "MSG")) +
  theme(legend.title=element_blank())+
  ggtitle("MSG vs PSG") + theme(plot.title = element_text(hjust = 0.5)) +
  #geom_text(aes(label =Geneid), color = "black", size = 3)
  
  geom_text(aes(label = ifelse(-log10(padj) > threshold, Geneid, "")),
            color = "black", size = 3, nudge_y = 4)  
vol_plot  

#Add the XP for each LOC
library(readxl)
Geneid_XP <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx") #gene annotation list 
Geneid_XP <- Geneid_XP %>%
  select('Geneid', 'Protein accession') %>%
  drop_na()

All_Plodia_genes_for_flybase_blastp <-Geneid_XP %>%
  drop_na()  %>%
  select(`Protein accession`)
write.table(All_Plodia_genes_for_flybase_blastp, file="All_Plodia_genes_for_flybase_prot_accessions.txt", quote=F, sep="\t",row.names=FALSE)


# input accessions in protein database https://www.ncbi.nlm.nih.gov/sites/batchentrez to get fasta with amino acid sequences 
flybase_results <- read_tsv("all_Plodia_proteins_for_flybase.results.txt", col_names=FALSE)
flybase_prot_to_Symbol <- read_tsv("dmel_unique_protein_isoforms_fb_2025_01.tsv", col_names=TRUE, comment="#")
flybase_gn_to_bpp <- read_tsv("fbgn_fbtr_fbpp_expanded_fb_2025_01.tsv", col_names=TRUE, comment="#") 
flybase_gn_summary <- read_tsv("automated_gene_summaries_fb_2025_01.tsv", col_names=TRUE, comment="#") 


flybase_prot_to_Symbol<- flybase_prot_to_Symbol %>%
  select(FBgn, FB_gene_symbol) %>%
  distinct(FBgn, .keep_all = TRUE) %>%
  left_join(flybase_gn_summary, by = c("FBgn" = "FBgn"))

flybase_gn_to_bpp <- flybase_gn_to_bpp %>%
  select(gene_ID, gene_fullname, polypeptide_ID) %>%
  distinct(polypeptide_ID, .keep_all = TRUE) %>%
  rename("gene_ID" = "FlyBase_FBgn") %>%
  drop_na()

flybase_results <- flybase_results %>%
  select(X1,X2,X11)

colnames(flybase_results)  <- c("Protein", "Flybase_prot", "Eval")

flybase_results <- flybase_results %>%
  left_join(Geneid_XP, by = c("Protein" = "Protein accession")) %>%
  relocate("Geneid", .after = "Protein") %>%
  group_by(Geneid) %>%
  slice_min(order_by = .data$Eval, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(flybase_gn_to_bpp, by = c("Flybase_prot" = "polypeptide_ID")) %>%
  left_join(flybase_prot_to_Symbol, by = c("FlyBase_FBgn" = "FBgn")) 

flybase_results_clean <- flybase_results %>%
  select(-Protein)

Plodia_gene_IDs <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx", col_names=TRUE) 
Plodia_gene_IDs <- Plodia_gene_IDs %>%
  select(Geneid, Symbol, Name) %>%
  rename("Name" = "NCBI_Name")  %>%
  distinct(Symbol, .keep_all = TRUE)


topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Geneid")) %>%
  relocate("NCBI_Name", .after = "Geneid")

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  relocate("Symbol", .after = "Geneid")



topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Geneid"))


write.table(topGenesMSG_PSG, file="DESeq2_Results_MSG_vs_PSG.tsv", quote=F, sep="\t",row.names=FALSE, na="")



# to get the normalized counts per tissue 

counts_ddsObj <- counts(ddsObj, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Geneid")


counts_ddsObj <- counts_ddsObj %>%
  mutate(MSG_avg = rowMeans(select(., ends_with("A")))) %>%
  mutate(PSG_avg = rowMeans(select(., ends_with("B"))))%>%
  mutate(SalG_avg = rowMeans(select(., ends_with("C"))))%>%
  mutate(Head_avg = rowMeans(select(., ends_with("D")) %>% select(where(is.numeric)))) %>%
  select("Geneid", ends_with("avg"), ends_with("A"), ends_with("B"), ends_with("C"), ends_with("D"),  everything()) 

counts_ddsObj <- counts_ddsObj %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Geneid")) %>%
  relocate("NCBI_Name", .after = "Geneid")

counts_ddsObj <- counts_ddsObj %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Geneid"))%>%
  relocate("MSG_avg", "PSG_avg", "SalG_avg", "Head_avg", "1A", "1B", "1C", "1D", "2A", "2B", "2C", "2D", "3A", "3B", "3C", "3D", "4A", "4B", "4C", "4D", .after = "Geneid") %>%
  relocate("NCBI_Name", .after = "Geneid") 

write.table(counts_ddsObj, file="All_Plodia_SG_DESeq2_normalized_counts.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# Now calculating TPM 

# Read the counts file in
seqdata <- read_tsv("Pi_SG_SalG.new.annotation.featurecounts.txt", comment="#")

# Add your curated gene symbols instead 
#first add Geneid for each old symbol and then add the new symbol
Symbol_to_Geneid <- read_excel("Plodia_gene_IDs.xlsx")

Symbol_to_Geneid <- Symbol_to_Geneid %>%
  select(Geneid,Symbol)%>%
  distinct(Symbol, .keep_all = TRUE)

seqdata <- seqdata %>%
  rename("Geneid" = "old_symbol") 

manual_annotations <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx")

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

seqdata <- seqdata %>%
  left_join(Symbol_to_Geneid, by = c("old_symbol" = "Symbol"))  %>%
  drop_na()


seqdata <- seqdata %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid"))  %>%
  select(-`old_symbol`) 



# Simplify to sample name for each counts column 

seqdata <- seqdata %>%
  rename_with(~ ifelse(
    str_starts(., "/"),
    str_extract(., "[^/]+(?=_pass2_mappedAligned\\.sortedByCoord\\.out\\.bam)"),
    .
  ))
# get gene lengths 
gene_lengths <- seqdata %>%
  dplyr::select(Geneid, Length)

# Convert gene_lengths to a named vector
gene_lengths <- setNames(gene_lengths$Length, gene_lengths$Geneid)

# Remove columns other than GeneID and sample counts
seqdata = seqdata %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

# Add the true NCBI gene id for each gene not the NCBI symbol by retrieving them from the NCBI spreadsheet but add "LOC" to Geneid number except for mitochondrial genes
Symbol_to_Geneid <- read_excel("Plodia_gene_IDs.xlsx")
Symbol_to_Geneid <- Symbol_to_Geneid %>%
  select(Geneid,Symbol)

# Transform raw data to matrix of counts:
countdata <- as.data.frame(seqdata)  %>%
  group_by(Geneid) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  column_to_rownames("Geneid") %>%
  as.matrix() 


counts_to_tpm <- function(counts, gene_lengths) {
  # Convert lengths from bp to kilobases
  gene_lengths_kb <- gene_lengths / 1000
  
  # Calculate Reads Per Kilobase (RPK)
  rpk <- counts / gene_lengths_kb
  
  # Compute scaling factor (per sample)
  scaling_factors <- colSums(rpk)
  
  # Calculate TPM
  tpm <- t( t(rpk) / scaling_factors ) * 1e6
  
  return(tpm)
}

tpm_matrix <- counts_to_tpm(countdata, gene_lengths)

# reformatting of the result table 
tpm_matrix <- as.data.frame(tpm_matrix) %>%
  rownames_to_column("Geneid")

tpm_matrix <- tpm_matrix %>%
  dplyr::select(-`1B_cat`, -`2B_cat`, -`3B_cat`, -`4B_cat`)

# To sum the technical replicates in the PSG 

tpm_matrix <- tpm_matrix %>%
  rowwise() %>%
  mutate(
    `1B` = sum(c_across(any_of(c("1B_first_run", "1B"))), na.rm = TRUE) / 2,
    `2B` = sum(c_across(any_of(c("2B_first_run", "2B"))), na.rm = TRUE) / 2,
    `3B` = sum(c_across(any_of(c("3B_first_run", "3B"))), na.rm = TRUE) / 2,
    `4B` = sum(c_across(any_of(c("4B_first_run", "4B"))), na.rm = TRUE) / 2
  ) %>%
  ungroup() %>%
  select(-any_of(c(
    "1B_first_run", "2B_first_run", "3B_first_run", "4B_first_run"
  ))) %>%
  select(Geneid, everything())

tpm_matrix <- tpm_matrix %>%
  mutate(
    MSG_avg = rowMeans(across(ends_with("A"), .names = "MSG_{.col}"), na.rm = TRUE),
    PSG_avg = rowMeans(across(ends_with("B"), .names = "PSG_{.col}"), na.rm = TRUE),
    SalG_avg = rowMeans(across(ends_with("C"), .names = "SalG_{.col}"), na.rm = TRUE),
    Head_avg = rowMeans(across(ends_with("D") & where(is.numeric)), na.rm = TRUE)
  ) %>%
  dplyr::select(Geneid, ends_with("avg"), ends_with("A"), ends_with("B"), ends_with("C"), ends_with("D"), everything())

# input accessions in protein database https://www.ncbi.nlm.nih.gov/sites/batchentrez to get fasta with amino acid sequences 
flybase_results <- read_tsv("all_Plodia_proteins_for_flybase.results.txt", col_names=FALSE)
flybase_prot_to_Symbol <- read_tsv("dmel_unique_protein_isoforms_fb_2025_01.tsv", col_names=TRUE, comment="#")
flybase_gn_to_bpp <- read_tsv("fbgn_fbtr_fbpp_expanded_fb_2025_01.tsv", col_names=TRUE, comment="#") 
flybase_gn_summary <- read_tsv("automated_gene_summaries_fb_2025_01.tsv", col_names=TRUE, comment="#") 


flybase_prot_to_Symbol<- flybase_prot_to_Symbol %>%
  select(FBgn, FB_gene_symbol) %>%
  distinct(FBgn, .keep_all = TRUE) %>%
  left_join(flybase_gn_summary, by = c("FBgn" = "FBgn"))

flybase_gn_to_bpp <- flybase_gn_to_bpp %>%
  select(gene_ID, gene_fullname, polypeptide_ID) %>%
  distinct(polypeptide_ID, .keep_all = TRUE) %>%
  rename("gene_ID" = "FlyBase_FBgn") %>%
  drop_na()

flybase_results <- flybase_results %>%
  select(X1,X2,X11)

colnames(flybase_results)  <- c("Protein", "Flybase_prot", "Eval")

flybase_results <- flybase_results %>%
  left_join(Geneid_XP, by = c("Protein" = "Protein accession")) %>%
  relocate("Geneid", .after = "Protein") %>%
  group_by(Geneid) %>%
  slice_min(order_by = .data$Eval, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(flybase_gn_to_bpp, by = c("Flybase_prot" = "polypeptide_ID")) %>%
  left_join(flybase_prot_to_Symbol, by = c("FlyBase_FBgn" = "FBgn")) 

flybase_results_clean <- flybase_results %>%
  select(-Protein)

Plodia_gene_IDs <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx", col_names=TRUE) 
Plodia_gene_IDs <- Plodia_gene_IDs %>%
  select(Geneid, Symbol, Name) %>%
  rename("Name" = "NCBI_Name")  %>%
  distinct(Symbol, .keep_all = TRUE)

tpm_matrix  <- tpm_matrix  %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Geneid")) %>%
  relocate("NCBI_Name", .after = "Geneid") %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Geneid"))%>%
  relocate("FB_gene_symbol", .after = "Symbol")%>%
  relocate("gene_fullname", .after = "FB_gene_symbol") %>%
  dplyr::select(-Flybase_prot, -Eval, -FlyBase_FBgn, -summary_text) %>%
  relocate("Geneid", "NCBI_Name", .after = "Symbol") %>%
  relocate("Symbol", .after = "Geneid") 

write.table(tpm_matrix, file="all_genes_TPM_normalized_counts_final.tsv", quote=F, sep="\t", row.names=FALSE, na="")

# now make a volcano plot Log2TPM vs Log2FC 

allGenesMSGvsPSG <- allGenesMSG_PSG_dds %>%
  dplyr::select(Geneid, log2FoldChange, padj)

# add a column to the table "all_genes_TPM_normalized_counts_final.tsv" you just generated that is the Max TPM value in MSG or PSG and save then read that file in

TPM_normalized_counts <- read_excel("all_genes_TPM_normalized_counts_final.xlsx")
TPM_normalized_counts <- TPM_normalized_counts %>%
  select(Symbol, Geneid, `Max(MSG,PSG)`) %>%
  left_join(allGenesMSGvsPSG, by = c("Geneid" = "Geneid")) 


# define significance by padj < 0.05 and TPM > 0.1 
TPM_normalized_counts <- TPM_normalized_counts %>%
  
  mutate(direction = case_when(
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))

TPM_normalized_counts_filtered_byTPM <- TPM_normalized_counts %>%
  filter(log2(`Max(MSG,PSG)`) > -2)

library(ggrepel)

vol_plot <- TPM_normalized_counts_filtered_byTPM %>%
  ggplot(aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +  
  scale_y_continuous( limits=c(-2, 20), expand=c(0,0)) +
  scale_x_continuous( limits=c(-12, 12), expand=c(0,0)) +
  geom_point(data = filter(TPM_normalized_counts_filtered_byTPM, direction == "ns"),
             aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction),
             show.legend = TRUE) +
  geom_point(data = filter(TPM_normalized_counts_filtered_byTPM, direction != "ns"),
             aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +
  
  scale_color_manual(
    values = c("down" = "#26b3ff", "up" = "#bb0c00", "ns" = "grey"),
    labels = c("down" = "PSG", "up" = "MSG", "ns" = "Non-Significant")
  ) +
  theme_classic() +
  theme(legend.title=element_blank())+
  ggtitle("MSG vs PSG") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = ifelse((log2(`Max(MSG,PSG)`) > 10 & direction != "ns") | ((log2FoldChange) < -4 & direction != "ns") | ((log2FoldChange) > 6 & direction != "ns"), Symbol, "")),
                  color = "black", size = 3, nudge_y = 0.5, max.overlaps = 15) +
  ylab("log2TPM") 



vol_plot 


# now annotate only specified genes
manual_annotations <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx")

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

annotation_list <- trimws(readLines("volplot_annotation_list.txt"))

TPM_normalized_counts_filtered_byTPM <- TPM_normalized_counts_filtered_byTPM %>%
  filter(log2(`Max(MSG,PSG)`) > -2)

TPM_normalized_counts_filtered_byTPM  <- TPM_normalized_counts_filtered_byTPM  %>%
  mutate(direction = case_when(
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))

vol_plot <- TPM_normalized_counts_filtered_byTPM %>%
  ggplot(aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +  
  scale_y_continuous( limits=c(-2, 20), expand=c(0,0)) +
  scale_x_continuous( limits=c(-9, 12), expand=c(0,0)) +
  geom_point(data = filter(TPM_normalized_counts_filtered_byTPM, direction == "ns"),
             aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction),
             show.legend = TRUE) +
  geom_point(data = filter(TPM_normalized_counts_filtered_byTPM, direction != "ns"),
             aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +
  
  scale_color_manual(
    values = c("down" = "#26b3ff", "up" = "#bb0c00", "ns" = "grey"),
    labels = c("down" = "PSG", "up" = "MSG", "ns" = "Non-Significant")
  ) +
  theme_classic() +
  theme(legend.title=element_blank())+
  ggtitle("MSG vs PSG") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(
    aes(label = ifelse(Symbol %in% annotation_list, Symbol, "")),
    color = "black",
    size = 3,
    nudge_y = 0.7,                  # Nudges label vertically
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.3,
    fill = alpha("white", 0.7),
    label.r = unit(0, "lines"),     # No rounded corners
    label.size = 0,                 # No label border
    point.padding = 0,          # Ensures label does not overlap point
    box.padding = 0             # Controls spacing around labels
    
  ) +
  ylab("log2TPM") 

vol_plot

write.table(TPM_normalized_counts_filtered_byTPM, file="volplot_data.tsv", quote=F, sep="\t", row.names=FALSE, na="")

#HEATMAP CODE 

##############Final Plot collapsed by tissue#####################

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

#before reading in add Max(MSG, PSG) column
normcounts <- read_excel("All_Plodia_SG_DESeq2_normalized_counts.xls") 

normcounts <- normcounts %>%
  select(Geneid, `Max(MSG,PSG)`) 


topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(normcounts, by = c("Geneid" = "Geneid"))


# Filter significant genes (adjust as needed)
sigGenes <- data.frame(Geneid = topGenesMSG_PSG$Geneid[
  topGenesMSG_PSG$padj <= 0.01 & abs(topGenesMSG_PSG$log2FoldChange) > 3 & (topGenesMSG_PSG$`Max(MSG,PSG)`) > 300
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

# add your manyal annotations

manual_annotations <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx")

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

Z <- Z %>%
  left_join(manual_annotations, by = c("Geneid" = "Geneid")) %>%
  left_join(flybasenames_unique, by = c("Geneid" = "Geneid")) %>%
  mutate(
    Symbol = ifelse(str_detect(Symbol, "LOC") & !is.na(FB_gene_symbol), FB_gene_symbol, Symbol)
  ) %>%
  select(-FB_gene_symbol, -Geneid) 

library(reshape2); library(ggdendro)

Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")

Z_df <- Z_df %>%
  group_by(Gene, Sample) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

#convert to wide matrix format for clustering
Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL


# Compute distances and clusters
distanceGene <- dist(Z_df_matrix)
distanceSample <- dist(t(Z_df_matrix))
clusterSample <- hclust(distanceSample, method = "average")
clusterGene <- hclust(distanceGene, method = "average")


#make dendogram for genes 
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
geneDendrogram <- ggplot(geneDendrogramData) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dendro()

Z_df$Sample  <- factor(Z_df$Sample , levels=c("MSG","PSG","SalG","Head"))
Z_df$Gene <- factor(Z_df$Gene, levels=clusterGene$labels[clusterGene$order])

#define your color palette 
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm
library(gridExtra); library(patchwork); library(cowplot)

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
  
  
write.table(Z_df_matrix, file="heatmap_raw_table_229genes.tsv", quote=F, sep="\t",row.names=TRUE, na="")


