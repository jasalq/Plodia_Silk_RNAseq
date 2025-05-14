# Plodia_Silk_RNAseq
[Description for project with the pre-print link]

## Pipeline Outline 
1. [RNAseq Quality Assessment](#assessment-of-rna-sequencing-quality-using-fastqc) 
2. [Adapter Trimming](#adapter-trimming-using-fastp)
3. [Read Mapping to Reference Genome](#star-mapping-of-rnaseq-data-to-the-<em>Plodia</em>-reference-genome-(GCF-027563975.2))
4. Read Counting
5. Differential Expression Analysis (DESeq2)
6. Gene Identification through Homology Searches
7. Visualization of Differential Expression Analysis Results in R
    
### Tools Used 
* Blast+ v.2.16.0+
* Fastp v.0.21.0	
* FastQC v.0.11.8
* JDK v.21.0.1
* NCBI Datasets
* RStudio v2024.04.2+764
* R/ cowplot
* R/ DESeq2
* R/ dplyr
* R/ ggdendro
* R/ ggplot2
* R/ ggrepel
* R/ gridExtra
* R/ pals
* R/ patchwork
* R/ readxl
* R/ reshape2
* R/ stringr
* R/ tidyverse
* STAR v.2.7.11b
* Subread v.2.0.8

## Assessment of RNA Sequencing Quality using FastQC
Quality of RNA sequencing was accessed using FastQC
```
#!/bin/sh
#SBATCH -J Plodia_silk_gland_rnaseq_fastqc_first_Run
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o Plodia_silk_gland_rnaseq_fastqc_first_run.out #redirecting stdout
#SBATCH -p debug #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 38 #amount of cores 
#SBATCH -t 1:00:00


echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

cd /scratch/martinlab/jasmine/30-1133207013/00_fastq

module load jdk/21.0.1
module load fastQC/0.11.8

mkdir /scratch/martinlab/jasmine/30-1133207013/00_fastq/fastqc_output_first_run

for i in *001_first_run.fastq.gz; do
        fastqc -f fastq -o /scratch/martinlab/jasmine/30-1133207013/00_fastq/fastqc_output_first_run $i;
        done

```
## Adapter Trimming using Fastp
Trim the adapters and polyG tails with Fastp
```
#!/bin/sh
#SBATCH -J plodia_SG_rnaseq_adap_trim
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o adap_trim_%A_%a.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 00:30:00
#SBATCH --array=0-19%9


echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

module load fastp

cd /scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed

PREFIX=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/Plodia_SG_SalG_Azenta_RNAseq_March2025/00_fastq

names=($(cat ${PREFIX}/samples))

echo ${names[${SLURM_ARRAY_TASK_ID}]} 

 fastp -L -Q --detect_adapter_for_pe --trim_poly_g \
        -i ${PREFIX}/${names[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz -I ${PREFIX}/${names[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz -o ${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R1_001.fastq.gz  -O ${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R2_001.fastq.gz
```
### Check that the adapters and polyG tails were removed using FastQC
Run FastQC again to check all adapters and polyG tails are gone

```
#!/bin/sh
#SBATCH -J Plodia_silk_gland_rnaseq_fastqc_second_run
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o Plodia_silk_gland_rnaseq_fastqc_second_run.out #redirecting stdout
#SBATCH -p debug #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 38 #amount of cores 
#SBATCH -t 1:00:00


echo "=========================================================="
echo "Running on node : $HOSTNAME"
date
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

cd /scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed

module load jdk/21.0.1
module load fastQC/0.11.8


names=($(cat samples_R1_R2))

echo ${names[${SLURM_ARRAY_TASK_ID}]} 

mkdir /scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed/fastqc_output

fastqc -f fastq -o /scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed/fastqc_output ${names[${SLURM_ARRAY_TASK_ID}]}_001.fastq.gz 

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```

After the job is done rename the output files to be recongnizeable by the sample name 
```
for i in Plodia_silk_gland_rnaseq_fastqc_second_run_*.out; do
	line6=$(sed -n '6p' ${i})
	mv ${i} "fastqc_${line6}.out";
	done 
```
## STAR Mapping of RNAseq data to the <em>Plodia</em> reference genome (GCF_027563975.2)

**Download <em>Plodia</em> RefSeq genome from NCBI using NCBI Datasets tool** 
```
mamba activate NCBI_datsets
datasets download genome accession GCF_027563975.2 
unzip ncbi_dataset
mv ncbi_dataset/data/GCF_027563975.2/GCF_027563975.2_ilPloInte3.2_genomic.fna ../../../
```
**Create a STAR genome index for the <em>Plodia</em> genome**  
Prior to running, to determine the correct value for --genomeSAindexNbases I used the formula provided where the genome is 291.3 Mb: min(14, log2(291,300,000) / 2 - 1).

```
#!/bin/sh
#SBATCH -J plodia_genome_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o plodia_genome_index.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -N 1 #amount of nodes 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "=========================================================="

PATH=$PATH:/CCAS/groups/martinlab/jasmine/software/STAR/source

cd /scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/plodia_genome_index
GENOME=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2_ilPloInte3.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2.gtf

mkdir $GENOME_DIR
STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir $GENOME_DIR --genomeFastaFiles $GENOME --sjdbGTFfile $ANNOTATION --sjdbOverhang 149 --genomeSAindexNbases 13
```

**Run the first pass of STAR mapping to the <em>Plodia</em> genome** 
```
#!/bin/sh
#SBATCH -J plodia_SG_rnaseq_star_pass_1
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o star_pass_1_%A_%a.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00
#SBATCH --array=0-19%10

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

PATH=$PATH:/CCAS/groups/martinlab/jasmine/software/STAR/source

ulimit -n 10000

PREFIX=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs
names=($(cat ${PREFIX}/samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]} 

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/plodia_genome_index_out
GENOME=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2_ilPloInte3.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2.gtf
RNAseq_FILES_PATH=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed
OUT_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/star_pass1

cd ${OUT_DIR}

STAR --runMode alignReads --runThreadN $CORES --genomeDir $GENOME_DIR --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${names[${SLURM_ARRAY_TASK_ID}]} --readFilesCommand zcat \
        --sjdbGTFfile $ANNOTATION --limitBAMsortRAM 60000000000 \
        --sjdbOverhang 149 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --sjdbGTFtagExonParentTranscript Parent \
        --readFilesIn ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R1_001.fastq.gz ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R2_001.fastq.gz 
      
echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```

After job is finished rename each log file by sample name
```
for i in star_pass_1_*.out; do
	line8=$(sed -n '8p' ${i})
	mv ${i} "star_pass_1_${line8}.out";
	done 
```

Due to the repeated sequencing of PSG samples there were two seperate sequencing files for each sample. To see if concatenating these samples prior to differential expression analyses was the proper way to do the analysis, I performed the same first round of STAR mapping with files concatenated by sample.

**Concatenate files from the same sample that were resequenced**
```
cat 1B*R1_001.fastq.gz > 1B_cat_adap_trim_R1_001.fastq.gz 
cat 1B*R2_001.fastq.gz > 1B_cat_adap_trim_R2_001.fastq.gz 
cat 2B*R1_001.fastq.gz > 2B_cat_adap_trim_R1_001.fastq.gz 
cat 2B*R2_001.fastq.gz > 2B_cat_adap_trim_R2_001.fastq.gz 
cat 3B*R1_001.fastq.gz > 3B_cat_adap_trim_R1_001.fastq.gz 
cat 3B*R2_001.fastq.gz > 3B_cat_adap_trim_R2_001.fastq.gz 
cat 4B*R1_001.fastq.gz > 4B_cat_adap_trim_R1_001.fastq.gz 
cat 4B*R2_001.fastq.gz > 4B_cat_adap_trim_R2_001.fastq.gz
```
**Re-run the first pass of STAR alignment to the <em>Plodia</em> genome using with files**
```

#!/bin/sh
#SBATCH -J plodia_SG_rnaseq_star_pass_1_cat
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o star_pass_1_cat_%A_%a.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00
#SBATCH --array=0-3%4

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

PATH=$PATH:/CCAS/groups/martinlab/jasmine/software/STAR/source

ulimit -n 10000

PREFIX=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs
names=($(cat ${PREFIX}/cat_samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]} 

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/plodia_genome_index_out
GENOME=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2_ilPloInte3.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2.gtf
RNAseq_FILES_PATH=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed
OUT_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/star_pass1

cd ${OUT_DIR}

STAR --runMode alignReads --runThreadN $CORES --genomeDir $GENOME_DIR --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${names[${SLURM_ARRAY_TASK_ID}]} --readFilesCommand zcat \
        --sjdbGTFfile $ANNOTATION --limitBAMsortRAM 60000000000 \
        --sjdbOverhang 149 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --sjdbGTFtagExonParentTranscript Parent \
        --outReadsUnmapped Fastx \
        --readFilesIn ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R1_001.fastq.gz ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R2_001.fastq.gz 
      
echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
To rename files by sample rather than array job ID 
```
for i in star_pass_1_cat*.out; do
	line8=$(sed -n '8p' ${i})
	mv ${i} "star_pass_1_cat${line8}.out";
	done
```
**Run the second pass of STAR mapping to the <em>Plodia</em> genome**
```
mkdir star_pass2

#!/bin/sh
#SBATCH -J plodia_SG_rnaseq_star_pass_2
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o star_pass_2_%A_%a.out #redirecting stdout
#SBATCH -p nano #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 0:30:00
#SBATCH --array=0-23%10

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

PATH=$PATH:/CCAS/groups/martinlab/jasmine/software/STAR/source

ulimit -n 10000

PREFIX=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs
names=($(cat ${PREFIX}/all_samples))
echo ${names[${SLURM_ARRAY_TASK_ID}]} 

CORES=40
GENOME_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/plodia_genome_index_out
GENOME=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2_ilPloInte3.2_genomic.fna
ANNOTATION=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2.gtf
RNAseq_FILES_PATH=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/trimmed
OUT_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/star_pass2
PASS1_DIR=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/star_pass1

STAR --runMode alignReads --runThreadN $CORES --genomeDir $GENOME_DIR --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${names[${SLURM_ARRAY_TASK_ID}]}_pass2_mapped --readFilesCommand zcat \
        --sjdbGTFfile $ANNOTATION --limitBAMsortRAM 60000000000 \
        --sjdbOverhang 149 --outSAMstrandField intronMotif --alignSoftClipAtReferenceEnds No --sjdbGTFtagExonParentTranscript Parent \
        --outReadsUnmapped Fastx \
        --quantMode GeneCounts \
        --sjdbFileChrStartEnd ${PASS1_DIR}/${names[${SLURM_ARRAY_TASK_ID}]}_STARgenome/sjdbList.out.tab \
        --readFilesIn ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R1_001.fastq.gz ${RNAseq_FILES_PATH}/${names[${SLURM_ARRAY_TASK_ID}]}_adap_trim_R2_001.fastq.gz 
      

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
````
After the job is done rename the log files by sample rather than array job ID
```
for i in star_pass_2*.out; do
	line8=$(sed -n '8p' ${i})
	mv ${i} "star_pass_2_${line8}.out";
	done 
```
## Read Counts with FeatureCounts
Install subread software which contains FeatureCounts

```
conda create -n subread_2.0.8 -c conda-forge mamba
conda activate subread_2.0.8
conda install -c conda-forge -c bioconda subread
```
Activate the subread enviroment and submit the counts job using a custom GTF annotation file that contains manual annotations for silk gland related genes.
```
#!/bin/sh
#SBATCH -J plodia_SG_rnaseq_counts
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o counts.out #redirecting stdout
#SBATCH -p defq #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 4:00:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

ANNOTATION=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/GCF_027563975.2_Ser1B-edited.gtf #use gff and are to see if the results are the same
OUT=Pi_SG_SalG.featurecounts.txt
BAM_FILES=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/star_runs/star_pass2
CORES=40

featureCounts -a $ANNOTATION -o $OUT  -T $CORES -p --primary -t exon -g gene_id ${BAM_FILES}/*.bam 2> FeatCounts.log

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
## Differential Expression Analysis with DESeq2 in RStudio
An R script with all of the code detailed below is [here](https://github.com/jasalq/Plodia_Silk_RNAseq/blob/afa0d7c60f025473a66a0297f3f168ada8b567db/R_Scripts/MSG_vs_PSG_DeSeq.R)
Some code inspiration was taken from Luca Livraghi's code from [Tendolkar et al. 2024](https://doi.org/10.7554/eLife.90846.3)

Set working directory
```
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/Jasmine/Pi_SG_SalG_RNAseq_results/Plodia_SG_Diff_Exp")
```
Install and load necessary packages
```
install.packages("tidyverse"); install.packages("stringr"); install.packages("dplyr"); install.packages("DESeq2")
library(tidyverse); library(stringr); library(dplyr); library(DESeq2); library(ggplot2)

```
```
Read the sample information into a data frame
sampleinfo <- read_tsv("all_samples2.txt")
sampleinfo

 Sample       Tissue
   <chr>        <chr> 
 1 1A           MSG   
 2 1B           PSG   
 3 1B_first_run PSG   
 4 1C           SalG  
 5 1D           Head  
 6 2A           MSG   
 7 2B           PSG   
 8 2B_first_run PSG   
 9 2C           SalG  
10 2D           Head  
```

Read the counts file in
```
seqdata <- read_tsv("Pi_SG_SalG.new.annotation.featurecounts.txt", comment="#")
```
Remove columns other than GeneID and sample counts
```
seqdata = seqdata %>%
  select(-Chr, -Start, -End, -Strand, -Length)
```
Simplify to sample name for each counts column 
```
 seqdata <- seqdata %>%
  rename_with(~ ifelse(. == "Sample_ID", ., str_extract(., "(?<=star_pass2_new_annotation/)[^/]+?_pass2")), -1) %>%
  rename_with(~ str_remove(., "_pass2$"), -1)
```
Adding identifiable names to previously annotated silk and AMP genes 
```
LOC_geneid <- read_tsv("GeneID_to_LOC.txt")

seqdata <- seqdata %>%
  left_join(LOC_geneid, by = c("Geneid" = "LOC")) %>%  # Join on Geneid and LOC
  mutate(
    Geneid = ifelse(!is.na(prot), prot, Geneid)  # Replace Geneid with prot if there's a match
  ) %>%
  select(-prot) 
```
Transform raw data to matrix of counts:
```
countdata <- as.data.frame(seqdata)  %>%
  column_to_rownames("Geneid") %>%
  as.matrix()
```
Remove the genes that are not/lowly expressed:
```
keep <- rowSums(countdata) > 0 
head(keep)
table(keep)

countdata <-countdata[keep, ]
dim(countdata)
head(countdata)
```
QC of counts; count distributions and visualize
```
summary(countdata)
boxplot(countdata, las=2)
```
Look at variance-mean relationship
```
plot(x=rowMeans(countdata), y=rowSds(countdata),
     main="sd vs mean",
     xlim=c(0,10000),
     ylim=c(0, 5000))
abline(lm(rowSds(countdata) ~ rowMeans(countdata)), col ="magenta")
```
### DESeq2 Analysis 
Define %out% function 
```
`%out%` <- function(x, y) !(x %in% y)
```
Remove the concatenated read samples from samples and count data
```
sampleinfo <- sampleinfo %>%
  filter(.[[1]] %out% c("1B_cat", "2B_cat", "3B_cat", "4B_cat"))

countdata <- countdata[, !colnames(countdata) %in% c("1B_cat", "2B_cat", "3B_cat", "4B_cat")]
```
Adding identifier to collapse by sequencing runs 
```
for(i in 1:nrow(sampleinfo)){
sampleinfo$seq_run[i] <- strsplit(sampleinfo$Sample[i], '_')[[1]][1]}
```
Define DESeq variables 
```
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sampleinfo,
  design = ~ Tissue
)
```
Collapse sequencing runs counts by sample and then generate ddsObj by running DESeq
```
ddscoll <- collapseReplicates(dds, dds$seq_run, renameCols=TRUE) #collapsing sequencing runs of PSG
ddsObj.raw <- DESeq(ddscoll)
ddsObj <- estimateSizeFactors(ddsObj.raw)
colData(ddsObj.raw)
colData(ddsObj)
ddsObj <- DESeq(ddsObj.raw)
res <- results(ddsObj, alpha=0.05) #adding p-value cutoff
res
```
Display contrast names
```
resultsNames(ddsObj)
```
Get right contrasts, in this case MSG vs PSG
```
Tissue_MSG_vs_PSG <- results(ddsObj, alpha = 0.05, contrast=c("Tissue","MSG","PSG"))
```
Doing some count plots
```
plotCounts(ddsObj, "FibH", intgroup = "Tissue")
plotCounts(ddsObj, "Ser1A", intgroup = "Tissue")
```
Display most diff expressed genes for this contrast with a padj < 0.05 and generate a tsv table
```
sum(Tissue_MSG_vs_PSG$padj < 0.05, na.rm = TRUE)

topGenesMSG_PSG <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>% 
  head(2898)

allGenesMSG_PSG_dds <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) 
write.table(allGenesMSG_PSG_dds, file="All_Plodia_genes_MSG_PSG_dds.tsv", quote=F, sep="\t", row.names = FALSE)
```

Replace some of the LOC numbers with geneIDs
```
LOC_geneid <- read_tsv("GeneID_to_LOC.txt")

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(LOC_geneid, by = c("Geneid" = "LOC")) %>%  # Join on Geneid and LOC
  mutate(
   Geneid = ifelse(!is.na(prot), prot, Geneid)) %>%  # Replace Geneid with prot if there's a match
  select(-prot) 
```

########################## Need to reformat stuff below still

















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


