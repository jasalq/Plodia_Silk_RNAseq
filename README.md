# Plodia_Silk_RNAseq
[Description for project with the pre-print link]

## Pipeline Outline 
1. [RNAseq Quality Assessment](#assessment-of-rna-sequencing-quality-using-fastqc) 
2. [Adapter Trimming](#adapter-trimming-using-fastp)
3. [Read Mapping to Reference Genome](#star-mapping-of-rnaseq-data-to-the-<em>Plodia</em>-reference-genome-(GCF-027563975.2))
4. [Read Counting](#read-counts-with-featurecounts)
5. [Differential Expression Analysis (DESeq2)](#differential-expression-analysis-with-deseq2-in-rstudio)
6. [Gene Identification through Flybase Homology Searches](#homology-searches-using-reciprocal-blastp-to-flybase)
7. [Final Visualization of Differential Expression Analysis Results in R](#final-visualization-of-differential-expression-analysis-results-in-r)
    
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
Generate preliminary Log2FoldChange vs Padj Volcano Plot
```
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  mutate(direction = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "ns"
  ))
# trying to force significance level to be right in box plot but this doesn't work because the box plot code takes ddsObj
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

write.table(topGenesMSG_PSG, file="MSG_PSG_new_anno_results.txt", quote=F, sep="\t")
```
**If you want to get a table with the normalized counts per tissue** 
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
```

Filter for comparisons 
```
counts_ddsObj <- counts_ddsObj[counts_ddsObj$Geneid %in% topGenesMSG_PSG$Geneid, ]
write.table(counts_ddsObj, file="MSG_vs_PSG_4tissue_Plodia_normalized_counts.tsv", quote=F, sep="\t", row.names=FALSE, na="")
```


## Homology Searches Using Reciprocal BLASTp to Flybase
*Add the XP (protein sequence) for each gene in the Plodia genome and make a list of all the <em>Plodia</em> proteins to blast*	 
[Note: Plodia_gene_IDs.txt is downloaded from the NCBI RefSeq genome annotation package accession no. GCF_027563975.2]
```
Geneid_XP <- read_tsv("Plodia_gene_IDs.txt")
Geneid_XP <- Geneid_XP %>%
  select('Symbol','Protein accession') %>%
  drop_na()

All_Plodia_genes_for_flybase_blastp <-Geneid_XP %>%
  drop_na()  %>%
  select(`Protein accession`)
write.table(All_Plodia_genes_for_flybase_blastp, file="All_Plodia_genes_for_flybase_prot_accessions.txt", quote=F, sep="\t")
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
flybase_prot_to_Symbol <- read_tsv("dmel_unique_protein_isoforms_fb_2025_01.tsv", col_names=TRUE, comment="#")
flybase_gn_to_bpp <- read_tsv("fbgn_fbtr_fbpp_expanded_fb_2025_01.tsv", col_names=TRUE, comment="#") 
flybase_gn_summary <- read_tsv("automated_gene_summaries_fb_2025_01.tsv", col_names=TRUE, comment="#") 
```
Do some reformatting of the tables to extract the pertinent information
```
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
```
**Finally add the flybase information for each Plodia gene in your DESeq2 results matrix and output a final table**
```  
topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(flybase_results_clean, by = c("Geneid" = "Symbol"))%>%
  relocate("FB_gene_symbol", .after = "Geneid")%>%
  relocate("gene_fullname", .after = "FB_gene_symbol") %>%
  left_join(Plodia_gene_IDs, by = c("Geneid" = "Symbol")) %>%
  relocate("NCBI_Name", .after = "Geneid")

write.table(topGenesMSG_PSG, file="topGenesMSG_PSG_FINAL_Seq_runs_comb_4tissue.tsv", quote=F, sep="\t",row.names=FALSE, na="")
```

I still have a lot of uncharacterized genes so I am trying to get info for uncharacterized genes from the protein name description column from the <em>Plodia</em> genome annotation information table 
```
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
```
## Final Visualization of Differential Expression Analysis Results in R
### Generating a Heatmap with the expression z-score for the top 230 genes in the MSG vs PSG contrast for all 4 tissues 
(based on the tutorial here https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/) 	
An R script with all of the code detailed below is [here](#R_Scripts/heatmap_code_cleaned.R)

Load libraries
```
library(readxl);library(DESeq2);library(reshape2);library(ggdendro);library(khroma);library(viridis);library(gridExtra);library(cowplot);library(pals)
```
Perform variance stabilizing transformation
```
vsd <- vst(ddsObj)
```
Convert to z-scores
```
Z <- t(scale(t(assay(vsd))))

# Convert to dataframe and keep gene IDs
Z <- as.data.frame(Z) %>%
  rownames_to_column("Geneid")
```
Add the Max average normalized counts PSG, MSG for filtering 
```
setwd("/Users/jasminealqassar/Library/CloudStorage/GoogleDrive-j.alqassar@gwmail.gwu.edu/.shortcut-targets-by-id/1Hi7WIp_ha7vnyQUclvt-AJRl6AF8y1hj/Martin Lab/Jasmine/Pi_SG_SalG_RNAseq_results/Plodia_SG_Diff_Exp")

normcounts <- read_excel("Plodia_SG_SalG_Diff_Exp_normalized_counts.xlsx", sheet = "MSG_vs_PSG_diff_exp_4tissue") 

normcounts <- normcounts %>%
  select(Geneid, `max(msg,psg)`)

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(normcounts, by = c("Geneid" = "Geneid"))
```

Filter significant genes (adjust as needed)
```
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

write.table(flybasenames_unique, file="flybasenames.tsv", quote=F, sep="\t", row.names=FALSE)

Z <- Z %>%
  left_join(flybasenames_unique, by = c("Geneid" = "Geneid")) %>%
  mutate(
    Geneid = ifelse(!is.na(FB_gene_symbol), FB_gene_symbol, Geneid)  # Replace Geneid with FB_symbol if match exists
  ) %>%
  select(-FB_gene_symbol) 
 ```
Melt to long format for ggplot2
```
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")
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

Re-factor samples for ggplot2 (only applied gene clustering in final figure)
```
#Z_df$Sample  <- factor(Z_df$Sample , levels=clusterSample$labels[clusterSample$order])
Z_df$Gene <- factor(Z_df$Gene, levels=clusterGene$labels[clusterGene$order])
```

Define your color palettes you might use
```
library(khroma)
vik <- color("vik")
vanimo <- color("vanimo")
batlowK <- color("batlowK")
tokyo <- color("tokyo")
oslo <- color("oslo")
berlin <- color("berlin")
library(viridis)
```

Construct the heatmap
```
heatmap <- ggplot(Z_df, aes(x=Sample, y=Gene, fill=Expression)) + geom_raster() + scale_fill_gradientn(colours =berlin(256)) + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y = element_blank(), axis.ticks.y=element_blank())
heatmap
grid.arrange(geneDendrogram, heatmap, ncol=1, heights=c(1,5))
```
Clean up the plot layout with cowplot 
```
# Remove duplicated legend
heatmap_clean <- heatmap + theme(legend.position = "none")

# Arrange plots
#top_row <- plot_grid(NULL, sampleDendrogram, ncol = 2, rel_widths = c(0.2, 1))
mid_row <- plot_grid(geneDendrogram, heatmap_clean, ncol = 2, rel_widths = c(0.2, 1))
legend <- get_legend(heatmap)

# Create heatmap
heatmap <- ggplot(Z_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_raster() +
  scale_fill_gradientn(colours = berlin(256)) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank()
  )
# Final layout
final_plot <- plot_grid(mid_row, ncol = 1, rel_heights = c(0.2, 1))
plot_grid(final_plot, legend, rel_widths = c(1, 0.12))
```

### Final Plot collapsed by tissue
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

# Convert to data frame
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

normcounts <- read_excel("Plodia_SG_SalG_Diff_Exp_normalized_counts.xlsx", sheet = "MSG_vs_PSG_diff_exp_4tissue") 

normcounts <- normcounts %>%
  select(Geneid, `max(msg,psg)`)

topGenesMSG_PSG <- topGenesMSG_PSG %>%
  left_join(normcounts, by = c("Geneid" = "Geneid"))
```

Filter significant genes (adjust as needed)
```
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
```
Melt to long format for ggplot2
```
Z_df <- melt(Z)
Z_df <- na.omit(Z_df)
colnames(Z_df) <- c("Gene", "Sample", "Expression")
```
Convert to wide matrix format for clustering
```
Z_df_matrix <- dcast(Z_df, Gene ~ Sample, value.var = "Expression")
rownames(Z_df_matrix) <- Z_df_matrix$Gene
Z_df_matrix$Gene <- NULL
```
Compute distances and clusters
```
distanceGene <- dist(Z_df_matrix)
distanceSample <- dist(t(Z_df_matrix))
clusterSample <- hclust(distanceSample, method = "average")
clusterGene <- hclust(distanceGene, method = "average")
```
Construct a dendogram for samples
```
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()
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
Re-factor samples for ggplot2 for ordering 
```
Z_df$Sample  <- factor(Z_df$Sample , levels=c("MSG","PSG","SalG","Head"))
Z_df$Gene <- factor(Z_df$Gene, levels=clusterGene$labels[clusterGene$order])
```
Define your color palette
```
library(pals)
ocean.balance <- ocean.balance
coolwarm <- coolwarm
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

# Now fix the heatmap and dendogram alignment issue
install.packages("patchwork")
library(patchwork)

# Remove legend from heatmap
heatmap_clean <- heatmap + theme(legend.position = "none")

# Combine dendrogram (left) + heatmap (right) side-by-side
main_panel <- geneDendrogram + heatmap_clean + plot_layout(ncol = 2, widths = c(1, 4))
main_panel <- geneDendrogram + heatmap_clean + 
  plot_layout(ncol = 2, widths = c(1, 6))  & 
  theme(plot.margin = margin(0, 0, 0, 0))

# Add legend to the right
final_plot <- main_panel + plot_spacer() + 
  plot_layout(ncol = 2, widths = c(5, 0.5))

legend <- get_legend(heatmap)
final_combined <- plot_grid(main_panel, legend, rel_widths = c(1, 0.12))

final_plot
```
