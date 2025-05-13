# Plodia_Silk_RNAseq
[Description for project with the pre-print link]

## Pipeline Outline 

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

