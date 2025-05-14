# Plodia Silk Gland Diff Exp Analysis (MSG vs PSG)
# Written by Jasmine D. Alqassar 2025 to use DESeq2 to analyze Plodia partitioned silk gland, salivary gland, and head RNAseq

# Set working directory
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/Jasmine/Pi_SG_SalG_RNAseq_results/Plodia_SG_Diff_Exp")

# Install necessary packages
#install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")

# Load necessary packages
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2)

# Read the sample information into a data frame
sampleinfo <- read_tsv("all_samples2.txt")
sampleinfo

# Read the counts file in
seqdata <- read_tsv("Pi_SG_SalG.new.annotation.featurecounts.txt", comment="#")

# Remove columns other than GeneID and sample counts
seqdata = seqdata %>%
  select(-Chr, -Start, -End, -Strand, -Length)

# Simplify to sample name for each counts column 
 seqdata <- seqdata %>%
  rename_with(~ ifelse(. == "Sample_ID", ., str_extract(., "(?<=star_pass2_new_annotation/)[^/]+?_pass2")), -1) %>%
  rename_with(~ str_remove(., "_pass2$"), -1)

#Adding identifiable names to previously annotated silk and AMP genes 
 
LOC_geneid <- read_tsv("GeneID_to_LOC.txt")

seqdata <- seqdata %>%
  left_join(LOC_geneid, by = c("Geneid" = "LOC")) %>%  # Join on Geneid and LOC
  mutate(
    Geneid = ifelse(!is.na(prot), prot, Geneid)  # Replace Geneid with prot if there's a match
  ) %>%
  select(-prot) 

# Transform raw data to matrix of counts:
countdata <- as.data.frame(seqdata)  %>%
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

# Define DESeq variables 

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sampleinfo,
  design = ~ Tissue
)
# Collapse sequencing runs counts by sample 
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

#doing some count plots
plotCounts(ddsObj, "FibH", intgroup = "Tissue")
plotCounts(ddsObj, "Ser1A", intgroup = "Tissue")

#Display most diff expressed genes for this contrast

sum(Tissue_MSG_vs_PSG$padj < 0.05, na.rm = TRUE)

topGenesMSG_PSG <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>% 
  head(2898)

allGenesMSG_PSG_dds <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) 
write.table(allGenesMSG_PSG_dds, file="All_Plodia_genes_MSG_PSG_dds.tsv", quote=F, sep="\t", row.names = FALSE)


# replace some of the LOC numbers with geneIDs
LOC_geneid <- read_tsv("GeneID_to_LOC.txt")

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(LOC_geneid, by = c("Geneid" = "LOC")) %>%  # Join on Geneid and LOC
  mutate(
   Geneid = ifelse(!is.na(prot), prot, Geneid)) %>%  # Replace Geneid with prot if there's a match
  select(-prot) 

#Volcano Plot 
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))
### trying to force significance level to be right in box plot but this doesn't work because the box plot code takes ddsObj
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "*",
    padj < 0.01 & log2FoldChange > 0 ~ "**",
    padj < 0.001 & log2FoldChange > 0 ~ "***",
    padj < 0.05 & log2FoldChange < 0 ~ "*",
    padj < 0.01 & log2FoldChange < 0 ~ "**",
    padj < 0.001 & log2FoldChange < 0 ~ "***",
    padj > 0.05 ~ "ns"
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

#write.table(topGenesMSG_PSG, file="MSG_PSG_new_anno_results.txt", quote=F, sep="\t")

#Add the XP for each LOC
Geneid_XP <- read_tsv("Plodia_gene_IDs.txt")
Geneid_XP <- Geneid_XP %>%
  select('Symbol','Protein accession') %>%
  drop_na()

All_Plodia_genes_for_flybase_blastp <-Geneid_XP %>%
  drop_na()  %>%
  select(`Protein accession`)
write.table(All_Plodia_genes_for_flybase_blastp, file="All_Plodia_genes_for_flybase_prot_accessions.txt", quote=F, sep="\t")


topGenesMSG_PSG_for_flybase <- topGenesMSG_PSG %>%
  left_join(Geneid_XP, by = c("Geneid" = "Symbol"))  # Join on Geneid and LOC

topGenesMSG_PSG_for_flybase <- topGenesMSG_PSG_for_flybase %>%
  drop_na() %>%
  select(`Protein accession`)# Join on Geneid and LOC

write.table(topGenesMSG_PSG_for_flybase, file="topGenesMSG_PSG_for_flybase_prot_accessions.txt", quote=F, sep="\t")

# edit the table to get rid of numbers in first column and header and input accessions in protein database https://www.ncbi.nlm.nih.gov/sites/batchentrez to get fasta with amino acid sequences 
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

  
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Symbol"))%>%
  relocate("FB_gene_symbol", .after = "Geneid")%>%
  relocate("gene_fullname", .after = "FB_gene_symbol") %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Symbol")) %>%
  relocate("NCBI_Name", .after = "Geneid")

write.table(topGenesMSG_PSG, file="topGenesMSG_PSG_FINAL_Seq_runs_comb_4tissue.tsv", quote=F, sep="\t",row.names=FALSE, na="")

# make a table with the differences to see what is happening 

df1 <- read_tsv("topGenesMSG_PSG_FINAL.tsv", col_names=TRUE) 
df2 <- read_tsv("topGenesMSG_PSG_FINAL_Seq_runs_comb.tsv", col_names=TRUE)

# Get rows unique to df1
only_df1 <- df1[!df1$Geneid %in% df2$Geneid, ]
only_df1$source <- "orig"

# Get rows unique to df2
only_df2 <- df2[!df2$Geneid %in% df1$Geneid, ]
only_df2$source <- "seq_runs_combined"

# Combine them
df3 <- bind_rows(only_df1, only_df2)

df3 <- df3 %>%
  relocate("source", .after = "Geneid")


write.table(df3, file="difference_between_combining_seqruns.tsv", quote=F, sep="\t", row.names=FALSE)

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
  left_join(LOC_geneid, by = c("Geneid" = "LOC")) %>%  # Join on Geneid and LOC
  mutate(
    Geneid = ifelse(!is.na(prot), prot, Geneid)  # Replace Geneid with prot if there's a match
  ) %>%
  select(-prot) 

counts_ddsObj <- counts_ddsObj %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Symbol"))%>%
  relocate("FB_gene_symbol", .after = "Geneid")%>%
  relocate("gene_fullname", .after = "FB_gene_symbol") %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Symbol")) %>%
  relocate("NCBI_Name", .after = "Geneid") %>%
  select(-Flybase_prot, -Eval, -FlyBase_FBgn, -summary_text)


write.table(counts_ddsObj, file="All_Plodia_SG_normalized_counts.tsv", quote=F, sep="\t", row.names=FALSE, na="")


#Filter for comparisons 

counts_ddsObj <- counts_ddsObj[counts_ddsObj$Geneid %in% topGenesMSG_PSG$Geneid, ]
write.table(counts_ddsObj, file="MSG_vs_PSG_4tissue_Plodia_normalized_counts.tsv", quote=F, sep="\t", row.names=FALSE, na="")

# generate original unfiltered counts_ddsObj and then run this 
SG_vs_SalG_manual_contrast <- read_tsv("topGenesSG_vs_SalG_manual_Contrast.tsv", col_names=TRUE, comment="#")
counts_ddsObj <- counts_ddsObj[counts_ddsObj$Geneid %in% SG_vs_SalG_manual_contrast$Geneid, ] #filtering to only include differentially expressed genes from SG vs SalG contrast
write.table(counts_ddsObj, file="topGenesSG_vs_SalG_manual_Contrast_normalized_counts.tsv", quote=F, sep="\t", row.names=FALSE, na="")


# trying to get info for uncharacterized genes 

topGenesMSGvsPSG <- read_excel("Plodia_SG_SalG_Diff_Exp_results_new_annotation.xlsx", sheet = "MSG_vs_PSG_4tissue") 

Geneid_XP <- read_tsv("Plodia_gene_IDs.txt")
Geneid_XP <- Geneid_XP %>%
  select(Symbol, `Protein name`)



uncharacterized_topGenesMSGvsPSG <- topGenesMSGvsPSG %>%
  filter(str_detect(NCBI_Name, "uncharacterized"))

uncharacterized_topGenesMSGvsPSG <- uncharacterized_topGenesMSGvsPSG %>%
left_join(Geneid_XP, by = c("Geneid" = "Symbol"))

uncharacterized_topGenesMSGvsPSG <- uncharacterized_topGenesMSGvsPSG %>%
  relocate("Protein name", .after = "NCBI_Name") %>%
  filter(!str_detect(`Protein name`, regex("uncharacterized", ignore_case = TRUE)))

uncharacterized_topGenesMSGvsPSG <- uncharacterized_topGenesMSGvsPSG %>%
  select(Geneid, `Protein name`)  %>%
  distinct(Geneid, .keep_all = TRUE)

topGenesMSGvsPSG <- read_excel("Plodia_SG_SalG_Diff_Exp_results_new_annotation.xlsx", sheet = "MSG_vs_PSG_4tissue") 

topGenesMSGvsPSG <- topGenesMSGvsPSG %>% 
  left_join(uncharacterized_topGenesMSGvsPSG, by = c("Geneid" = "Geneid")) %>%  # Join on Geneid and LOC
  mutate(
    NCBI_Name = ifelse(!is.na(`Protein name`), `Protein name`, NCBI_Name)  # Replace Geneid with prot if there's a match
  ) %>%
  select(-`Protein name`) 

write.table(topGenesMSGvsPSG, file="topGenesMSG_PSG_FINAL_Seq_runs_comb_4tissue.tsv", quote=F, sep="\t", row.names=FALSE, na="")
