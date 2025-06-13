setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/Jasmine/Pi_SG_SalG_RNAseq_results/Plodia_SG_Diff_Exp")


library(dplyr)
library(tidyverse)

# Read the counts file in
seqdata <- read_tsv("Pi_SG_SalG.new.annotation.featurecounts.txt", comment="#")

# get gene lengths 
gene_lengths <- seqdata %>%
  dplyr::select(Geneid, Length)

# Convert gene_lengths to a named vector
gene_lengths <- setNames(gene_lengths$Length, gene_lengths$Geneid)

# Remove columns other than GeneID and sample counts
seqdata = seqdata %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

# Simplify to sample name for each counts column 
seqdata <- seqdata %>%
  rename_with(~ ifelse(. == "Sample_ID", ., str_extract(., "(?<=star_pass2_new_annotation/)[^/]+?_pass2")), -1) %>%
  rename_with(~ str_remove(., "_pass2$"), -1)


# Transform raw data to matrix of counts:
countdata <- as.data.frame(seqdata)  %>%
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

# if you want to sum the technical replicates in the PSG 
tpm_matrix <- tpm_matrix %>%
  mutate(`1Bcoll` = rowSums(across(c(`1B_first_run`, `1B`)))) %>%
  mutate(`2Bcoll` = rowSums(across(c(`2B_first_run`, `2B`)))) %>%
  mutate(`3Bcoll` = rowSums(across(c(`3B_first_run`, `3B`)))) %>%
  mutate(`4Bcoll` = rowSums(across(c(`4B_first_run`, `4B`)))) %>%
  dplyr::select(-`1B`, -`2B`, -`3B`, -`4B`, -`1B_first_run`, -`2B_first_run`, -`3B_first_run`, -`4B_first_run` ) %>%
  dplyr::rename(
    `1B` = `1Bcoll`,
    `2B` = `2Bcoll`,
    `3B` = `3Bcoll`,
    `4B` = `4Bcoll`)  %>%
  select(order(colnames(tpm_matrix))) %>%
  dplyr::select(Geneid, everything())


LOC_geneid <- read_tsv("GeneID_to_LOC.txt")

tpm_matrix <- tpm_matrix %>%
  left_join(LOC_geneid, by = c("Geneid" = "LOC")) %>%  # Join on Geneid and LOC
  mutate(
    Geneid = ifelse(!is.na(prot), prot, Geneid)  # Replace Geneid with prot if there's a match
  ) %>%
  dplyr::select(-prot) 

tpm_matrix <- tpm_matrix %>%
  mutate(
    MSG_avg = rowMeans(across(ends_with("A"), .names = "MSG_{.col}"), na.rm = TRUE),
    PSG_avg = rowMeans(across(ends_with("B"), .names = "PSG_{.col}"), na.rm = TRUE),
    SalG_avg = rowMeans(across(ends_with("C"), .names = "SalG_{.col}"), na.rm = TRUE),
    Head_avg = rowMeans(across(ends_with("D") & where(is.numeric)), na.rm = TRUE)
  ) %>%
  dplyr::select(Geneid, ends_with("avg"), ends_with("A"), ends_with("B"), ends_with("C"), ends_with("D"), everything())

# edit the table to get rid of numbers in first column and header and input accessions in protein database https://www.ncbi.nlm.nih.gov/sites/batchentrez to get fasta with amino acid sequences 
flybase_results <- read_tsv("all_Plodia_proteins_for_flybase.results.txt", col_names=FALSE)
flybase_prot_to_Symbol <- read_tsv("dmel_unique_protein_isoforms_fb_2025_01.tsv", col_names=TRUE, comment="#")
flybase_gn_to_bpp <- read_tsv("fbgn_fbtr_fbpp_expanded_fb_2025_01.tsv", col_names=TRUE, comment="#") 
flybase_gn_summary <- read_tsv("automated_gene_summaries_fb_2025_01.tsv", col_names=TRUE, comment="#") 


flybase_prot_to_Symbol<- flybase_prot_to_Symbol %>%
  dplyr::select(FBgn, FB_gene_symbol) %>%
  distinct(FBgn, .keep_all = TRUE) %>%
  left_join(flybase_gn_summary, by = c("FBgn" = "FBgn"))

flybase_gn_to_bpp <- flybase_gn_to_bpp %>%
  dplyr::select(gene_ID, gene_fullname, polypeptide_ID) %>%
  distinct(polypeptide_ID, .keep_all = TRUE) %>%
  rename("gene_ID" = "FlyBase_FBgn") %>%
  drop_na()

flybase_results <- flybase_results %>%
  dplyr::select(X1,X2,X11)

colnames(flybase_results)  <- c("Protein", "Flybase_prot", "Eval")

Geneid_XP <- read_tsv("Plodia_gene_IDs.txt")
Geneid_XP <- Geneid_XP %>%
  dplyr::select('Symbol','Protein accession') %>%
  drop_na()


flybase_results <- flybase_results %>%
  left_join(Geneid_XP, by = c("Protein" = "Protein accession")) %>%
  relocate("Symbol", .after = "Protein") %>%
  group_by(Symbol) %>%
  slice_min(order_by = .data$Eval, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(flybase_gn_to_bpp, by = c("Flybase_prot" = "polypeptide_ID")) %>%
  left_join(flybase_prot_to_Symbol, by = c("FlyBase_FBgn" = "FBgn")) 

flybase_results_clean <- flybase_results %>%
  dplyr::select(-Protein)

Plodia_gene_IDs <- read_tsv("Plodia_gene_IDs.tsv", col_names=TRUE, comment="#") 
Plodia_gene_IDs <- Plodia_gene_IDs %>%
  dplyr::select(Symbol, Name) %>%
  rename("Name" = "NCBI_Name")  %>%
  distinct(Symbol, .keep_all = TRUE)

tpm_matrix  <- tpm_matrix  %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Symbol"))%>%
  relocate("FB_gene_symbol", .after = "Geneid")%>%
  relocate("gene_fullname", .after = "FB_gene_symbol") %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Symbol")) %>%
  relocate("NCBI_Name", .after = "Geneid") %>%
  dplyr::select(-Flybase_prot, -Eval, -FlyBase_FBgn, -summary_text)


write.table(tpm_matrix, file="all_genes_TPM_normalized_counts_PSG_coll.tsv", quote=F, sep="\t", row.names=FALSE, na="")


threshold <- 6 #assigning an arbitrary cutoff to label all points above this -log10(padj) value because the plot gets too crowded

library(readxl)
TPM_normalized_counts <- read_excel("all_genes_TPM_normalized_counts_PSG_coll.xlsx") 

topGenesMSGvsPSG <- read_tsv("topGenesMSG_PSG_FINAL_Seq_runs_comb_4tissue.tsv", col_names=TRUE, comment="#")

topGenesMSGvsPSG <- topGenesMSGvsPSG %>%
  dplyr::select(Geneid, log2FoldChange)

TPM_normalized_counts <- TPM_normalized_counts %>%
  dplyr::select(Geneid, `Max(MSG,PSG)`) %>%
  left_join(topGenesMSGvsPSG, by = c("Geneid" = "Geneid")) %>%
  drop_na()

# add flybase names for LOCs
topGenesMSGvsPSG <- read_tsv("topGenesMSG_PSG_FINAL_Seq_runs_comb_4tissue.tsv", col_names=TRUE, comment="#")

flybasenames <- topGenesMSGvsPSG %>%
  select(Geneid, FB_gene_symbol) %>%
  drop_na() %>%
  add_count(Geneid) %>%
  filter(n == 1) %>%
  select(-n)

flybasenames_unique <- flybasenames %>%
  group_by(FB_gene_symbol) %>%
  filter(n() == 1) %>%     # Keep only those with a single occurrence
  ungroup()


TPM_normalized_counts <- TPM_normalized_counts%>%
  left_join(flybasenames_unique, by = c("Geneid" = "Geneid")) %>%
  mutate(
    Geneid = ifelse(!is.na(FB_gene_symbol), FB_gene_symbol, Geneid)  # Replace Geneid with FB_symbol if match exists
  ) %>%
  select(-FB_gene_symbol) 

# define significance by padj < 0.05 and TPM > 0.1 
TPM_normalized_counts <- TPM_normalized_counts %>%
  mutate(direction = case_when(
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))






log2foldchange_threshold <- 5
# make a volcano plot 
library(ggrepel)
library(ggplot2)
vol_plot <- TPM_normalized_counts %>%
  ggplot(aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +  
  geom_point() +  
  scale_color_manual(values = c("down" = "#26b3ff", "ns" = "grey", "up" = "#bb0c00"),  
                     labels = c("down" = "PSG", "ns" = "Not significant", "up" = "MSG")) +
  theme(legend.title=element_blank()) +
  ggtitle("MSG vs PSG") + theme(plot.title = element_text(hjust = 0.5)) +
  
  geom_text_repel(aes(label = ifelse(log2(`Max(MSG,PSG)`) > threshold | abs(log2FoldChange) > log2foldchange_threshold, Geneid, "")),
            color = "black", size = 3, nudge_y = 0.5, max.overlaps = 15)  +
  
  ylab("log2TPM")

vol_plot  

# Now work on final plot with all genes plotted

library(readxl)
TPM_normalized_counts <- read_excel("all_genes_TPM_normalized_counts_PSG_coll.xlsx") 

allGenesMSGvsPSG <- read_tsv("All_Plodia_genes_MSG_PSG_dds.tsv", col_names=TRUE, comment="#")
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
  relocate("Symbol", .after = "Protein") %>%
  group_by(Symbol) %>%
  slice_min(order_by = .data$Eval, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(flybase_gn_to_bpp, by = c("Flybase_prot" = "polypeptide_ID")) %>%
  left_join(flybase_prot_to_Symbol, by = c("FlyBase_FBgn" = "FBgn")) 

flybase_results_clean <- flybase_results %>%
  select(-Protein)
Plodia_gene_IDs <- read_tsv("Plodia_gene_IDs.tsv", col_names=TRUE, comment="#") 
Plodia_gene_IDs <- Plodia_gene_IDs %>%
  select(Symbol, Name) %>%
  rename("Name" = "NCBI_Name")  %>%
  distinct(Symbol, .keep_all = TRUE)


allGenesMSGvsPSG <- allGenesMSGvsPSG %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Symbol"))%>%
  relocate("FB_gene_symbol", .after = "Geneid")%>%
  relocate("gene_fullname", .after = "FB_gene_symbol") %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Symbol")) %>%
  relocate("NCBI_Name", .after = "Geneid")

allGenesMSGvsPSG <- allGenesMSGvsPSG %>%
  dplyr::select(Geneid, log2FoldChange, padj)

TPM_normalized_counts <- TPM_normalized_counts %>%
  dplyr::select(Geneid, `Max(MSG,PSG)`) %>%
  left_join(allGenesMSGvsPSG, by = c("Geneid" = "Geneid")) %>%
  drop_na()

# define significance by padj < 0.05 and TPM > 0.1 
TPM_normalized_counts <- TPM_normalized_counts %>%
  mutate(direction = case_when(
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2(`Max(MSG,PSG)`) > 0.1 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))




threshold <- 1000
log2foldchange_threshold <- 5
# make a volcano plot 
library(ggrepel)
library(ggplot2)

vol_plot <- TPM_normalized_counts %>%
  ggplot(aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +  
  geom_point(data = filter(TPM_normalized_counts, direction == "ns"),
             aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction),
             show.legend = TRUE) +
  geom_point(data = filter(TPM_normalized_counts, direction != "ns"),
             aes(x = log2FoldChange, y = log2(`Max(MSG,PSG)`), color = direction)) +
  
  scale_color_manual(
    values = c("down" = "#26b3ff", "up" = "#bb0c00", "ns" = "grey"),
    labels = c("down" = "PSG", "up" = "MSG", "ns" = "Non-Significant")
  ) +
  
  theme_classic() +
  theme(legend.title=element_blank())+
  ggtitle("MSG vs PSG") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = ifelse((log2(`Max(MSG,PSG)`) > 10 & direction != "ns") | ((log2FoldChange) < -5 & direction != "ns") | ((log2FoldChange) > 7 & direction != "ns"), Geneid, "")),
                  color = "black", size = 3, nudge_y = 0.5, max.overlaps = 15) +
  ylab("log2TPM")
vol_plot 



write.table(TPM_normalized_counts, file="raw_TPM_foldchange_volplot_MSG_vs_PSG.xlsx", quote=F, sep="\t", row.names=FALSE, na="")

# FINAL PLOT Now trim the volplot so nothing below a log2TPM of -2 is plotted 

TPM_normalized_counts_filtered_byTPM <- TPM_normalized_counts %>%
filter(log2(`Max(MSG,PSG)`) > -2)


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
  geom_text_repel(aes(label = ifelse((log2(`Max(MSG,PSG)`) > 10 & direction != "ns") | ((log2FoldChange) < -5 & direction != "ns") | ((log2FoldChange) > 7 & direction != "ns"), Geneid, "")),
                  color = "black", size = 3, nudge_y = 0.5, max.overlaps = 15) +
  ylab("log2TPM") 



vol_plot 

write.table(TPM_normalized_counts_filtered_byTPM, file="final_filtered_June10_volplot_MSG_vs_PSG.xlsx", quote=F, sep="\t", row.names=FALSE, na="")

