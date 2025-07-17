# Plodia_Silk_RNAseq
<body>  
This repository contains all code associated with the manuscript Alqassar et al. 2025 <a href="https://doi.org/10.1101/2025.07.11.664249"> Regionalization of gene expression and cell types in the silk gland of the pantry moth <em>Plodia interpunctella </em> </a>. The code found in this repository processes, maps, and performs differential expression analyses on RNA sequencing data generated during the study where we aimed to understand differences in gene expression between different sections of silk gland tissue.		
<body/>	
<br><br>	
<strong>Contact info</strong>: Jasmine Alqassar (j.alqassar@gwu.edu) and Arnaud Martin (arnaud@gwu.edu)
<br><br>
<strong>Data Availability</strong>: All sequencing data generated during this study has been deposited in the NCBI Sequence Read Archive under BioProject PRJNA1241317. 

## Pipeline Outline 
1. [RNAseq Quality Assessment](#assessment-of-rna-sequencing-quality-using-fastqc) 
2. [Adapter Trimming](#adapter-trimming-using-fastp)
3. [Read Mapping to the Reference Genome](#star-mapping-of-rnaseq-data-to-the-plodia-reference-genome-gcf_0275639752)
5. [Read Counting](#read-counts-with-featurecounts)
6. [Differential Expression Analysis (DESeq2)](#differential-expression-analysis-with-deseq2-in-rstudio)
7. [Gene Identification through FlyBase Homology Searches](#homology-searches-using-reciprocal-blastp-to-flybase)
8. [TPM Normalization](#tpm-normalization)
9. [Heatmap Visualization of Differential Expression Analysis Results in R](#heatmap-visualization-of-differential-expression-analysis-results-in-r)
    
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
Add your curated gene symbols instead 
```
#First add Geneid for each old symbol and then add the new symbol

Symbol_to_Geneid <- read_excel("Plodia_gene_IDs.xlsx")

Symbol_to_Geneid <- Symbol_to_Geneid %>%
  select(Geneid,Symbol)%>%
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
```
Transform raw data to matrix of counts:
```
countdata <- seqdata  %>%
  group_by(Geneid) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
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
Define DESeq variables and design
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
Display most diff expressed genes for this contrast with a padj < 0.05 and generate a tsv table
```
sum(Tissue_MSG_vs_PSG$padj < 0.05, na.rm = TRUE)
topGenesMSG_PSG <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj) %>% 
  head(2914)

allGenesMSG_PSG_dds <- as.data.frame(Tissue_MSG_vs_PSG) %>%
  rownames_to_column("Geneid") %>% 
  arrange(padj)
```
A quick scatter plot
```
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))


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
```

## Homology Searches Using Reciprocal BLASTp to FlyBase
*Add the XP (protein sequence) for each gene in the Plodia genome and make a list of all the <em>Plodia</em> proteins to blast*	 
[Note: Plodia_gene_IDs.txt is downloaded from the NCBI RefSeq genome annotation package accession no. GCF_027563975.2]
```
#Add the XP for each LOC
Geneid_XP <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx") #gene annotation list 
Geneid_XP <- Geneid_XP %>%
  select('Geneid', 'Protein accession') %>%
  drop_na()

All_Plodia_genes_for_flybase_blastp <-Geneid_XP %>%
  drop_na()  %>%
  select(`Protein accession`)
write.table(All_Plodia_genes_for_flybase_blastp, file="All_Plodia_genes_for_flybase_prot_accessions.txt", quote=F, sep="\t",row.names=FALSE)
```
Visit [NCBI Batch Entrez](#https://www.ncbi.nlm.nih.gov/sites/batchentrez) and select protein database and upload the All_Plodia_genes_for_flybase_prot_accessions.txt you generated to retrieve a FASTA file with all the amino acid sequences for each of the proteins in the list you just generated

**Now run BLASTp on the HPC** 
```
#!/bin/sh
#SBATCH -J flybase_blastp_MSG_PSG
#SBATCH --mail-type=END
#SBATCH --mail-user=j.alqassar@gwu.edu
#SBATCH -o flybase_blastp_MSG_PSG.out #redirecting stdout
#SBATCH -p tiny;short #queue 
#SBATCH -n 40 #amount of cores 
#SBATCH -t 04:00:00

echo "=========================================================="
echo "Running on node : $HOSTNAME"
echo "Current directory : $PWD"
echo "Current job ID : $SLURM_JOB_ID"
echo "Job Started:"
date
echo "=========================================================="

module load blast+/2.16.0+

DB=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/get_flybase_names/db
QUERY=/scratch/martinlab/jasmine/plodia_silk_gland_diff_exp/get_flybase_names/All_Plodia_genes_for_flybase_prot.fa
PREFIX=all_Plodia_proteins_for_flybase

blastp -db ${DB} \
-query ${QUERY} \
-out ${PREFIX}.out \
-outfmt 6 \
-evalue 1e-5 \
-num_threads 40

echo "=========================================================="
echo "Job Finished  $SLURM_JOB_ID:"
date
echo "=========================================================="
```
**Now back to RStudio to add the results to your data matrix**	 

Read in the results
```
flybase_results <- read_tsv("all_Plodia_proteins_for_flybase.results.txt", col_names=FALSE)
```
Download and read in the following tables from the latest [FlyBase genome release](#https://flybase.org/)
```
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

```
### Now to add the normalized counts per tissue from DESeq2
```
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
```
## TPM Normalization
```
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/Jasmine/Pi_SG_SalG_RNAseq_results/Plodia_SG_Diff_Exp")

library(dplyr)
library(tidyverse)
```
Read the counts file in
```
seqdata <- read_tsv("Pi_SG_SalG.new.annotation.featurecounts.txt", comment="#")
```
Add your curated gene symbols instead of old NCBI symbols 
```
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

```
Simplify to sample name for each counts column 
```
seqdata <- seqdata %>%
  rename_with(~ ifelse(
    str_starts(., "/"),
    str_extract(., "[^/]+(?=_pass2_mappedAligned\\.sortedByCoord\\.out\\.bam)"),
    .
  ))
```
Get gene lengths 
```
gene_lengths <- seqdata %>%
  dplyr::select(Geneid, Length)
```
Convert gene_lengths to a named vector
```
gene_lengths <- setNames(gene_lengths$Length, gene_lengths$Geneid)
```
Remove columns other than GeneID and sample counts
```
seqdata = seqdata %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)
```
Add the true NCBI gene id for each gene not the NCBI symbol by retrieving them from the NCBI spreadsheet but add "LOC" to Geneid number except for mitochondrial genes
```
Symbol_to_Geneid <- read_excel("Plodia_gene_IDs.xlsx")
Symbol_to_Geneid <- Symbol_to_Geneid %>%
  select(Geneid,Symbol)
```
Transform raw data to matrix of counts:
```
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
```
Reformatting of the result table 
```
tpm_matrix <- as.data.frame(tpm_matrix) %>%
  rownames_to_column("Geneid")

tpm_matrix <- tpm_matrix %>%
  dplyr::select(-`1B_cat`, -`2B_cat`, -`3B_cat`, -`4B_cat`)
```
If you want to sum the technical replicates in the PSG 
```
# To sum the technical replicates for the PSG 
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
```
Now add your FlyBase matches to the TPM data 
```
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

```
### Now make a scatter plot Log2TPM vs Log2FC 
```
allGenesMSGvsPSG <- allGenesMSG_PSG_dds %>%
  dplyr::select(Geneid, log2FoldChange, padj)
```
Add a column to the table "all_genes_TPM_normalized_counts_final.tsv" you just generated that is the Max TPM value in MSG or PSG and save as excel sheet then read that file in
```
TPM_normalized_counts <- read_excel("all_genes_TPM_normalized_counts_final.xlsx")
TPM_normalized_counts <- TPM_normalized_counts %>%
  select(Symbol, Geneid, `"Max(MSG,PSG)"`) %>%
  rename(`"Max(MSG,PSG)"` = "Max(MSG,PSG)")  %>%
  left_join(allGenesMSGvsPSG, by = c("Geneid" = "Geneid")) %>%
  drop_na()
```

Define significance by padj < 0.05 and TPM > 0.1 
```
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
```

Now annotate only specified genes
```
manual_annotations <- read_excel("ilPloInte3.2_manually_curated_gene_table_final.xlsx")

manual_annotations <- manual_annotations %>%
  select(`Geneid`, Symbol)

annotation_list <- trimws(readLines("volplot_annotation_list.txt")) # this is a file I generated that is one column with the symbol names I wanted annotated

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
```

## Heatmap Visualization of Differential Expression Analysis Results in R 
(based on the tutorial here https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/) 	

Load libraries
```
library(readxl);library(DESeq2);library(reshape2);library(ggdendro);library(khroma);library(viridis);library(gridExtra);library(cowplot);library(pals)
```
Perform variance stabilizing transformation
```
vsd <- vst(ddsObj)
vst_mat <- assay(vsd)
```
Get tissue assignment for each sample
```
tissue_info <- colData(ddsObj)$Tissue
names(tissue_info) <- colnames(vst_mat)  # label with sample names
```
Create a list mapping each tissue to the sample columns
```
tissue_groups <- split(names(tissue_info), tissue_info)
```
Collapse vst matrix by computing rowMeans per tissue
```
collapsed_vst <- sapply(tissue_groups, function(samples) {
  rowMeans(vst_mat[, samples, drop = FALSE], na.rm = TRUE)
})
```
Convert to data frame
```
collapsed_vst_df <- as.data.frame(collapsed_vst) %>%
  rownames_to_column("Geneid")
```

Set rownames and remove Geneid column before scaling
```
mat <- collapsed_vst_df
rownames(mat) <- mat$Geneid
mat$Geneid <- NULL
```
Z-score across tissues (i.e., scale by row)
```
Z <- t(scale(t(as.matrix(mat))))

# Convert to dataframe and keep gene IDs
Z <- as.data.frame(Z) %>%
  rownames_to_column("Geneid")
```
Add the Max average normalized counts PSG, MSG for filtering 
```
normcounts <- read_excel("All_Plodia_SG_DESeq2_normalized_counts.xls") 

normcounts <- normcounts %>%
  select(Geneid, `Max(MSG,PSG)`) 

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(normcounts, by = c("Geneid" = "Geneid"))
```

Filter significant genes (adjust as needed)
```
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

# add your manual annotations
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
```

Melt to long format for ggplot2
```
library(reshape2); library(ggdendro)
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")

Z_df <- Z_df %>%
  group_by(Gene, Sample) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

```
Convert to wide matrix format for clustering
```
Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL
```
Compute distances and gene and sample clusters
```
distanceGene <- dist(Z_df_matrix)
distanceSample <- dist(t(Z_df_matrix))
clusterSample <- hclust(distanceSample, method = "average")
clusterGene <- hclust(distanceGene, method = "average")
```

Construct a dendogram for samples (commented out because we didn't end up clustering by sample in the final figure
```
#sampleModel <- as.dendrogram(clusterSample)
#sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
#sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()
```
Make dendogram for genes 
```
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
geneDendrogram <- ggplot(geneDendrogramData) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dendro()
```

Re-factor samples for ggplot2 
```
Z_df$Sample  <- factor(Z_df$Sample , levels=c("MSG","PSG","SalG","Head"))
Z_df$Gene <- factor(Z_df$Gene, levels=clusterGene$labels[clusterGene$order])

```

Define your color palettes you might use and load libraries
```
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm
library(gridExtra); library(patchwork); library(cowplot)
```
Construct the heatmap
```
heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =coolwarm(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))
```
Clean up the plot layout with cowplot 
```
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
```

