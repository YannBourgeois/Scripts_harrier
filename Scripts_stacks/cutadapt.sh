#!/bin/bash
#SBATCH --mem=12GB
#SBATCH --nodes=1
#SBATCH --job-name=stacks_job
#SBATCH --array=1-61         # 61 individuals 
#SBATCH --time=06:00:00
#SBATCH --output=stacks_job.%J.out
#SBATCH --error=stacks_job.%J.err
## ***Commands starts here***

config=list_for_cutadapt.txt

PATH_TO_SAMPLES=/mnt/lustre/yb24/Harrier_analyses/RE_processed/
PATH_TO_OUTPUT_SAMPLES=/mnt/lustre/yb24/Harrier_analyses/RE_processed/Second_reads_trimmed/
####This PATH can be changed to what works on your own cluster
ID_NUMBER=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
SAMPLE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
M_VALUE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

cutadapt -m 140 -l 140 -o ${PATH_TO_OUTPUT_SAMPLES}/${SAMPLE}.1.fq.gz -p  ${PATH_TO_OUTPUT_SAMPLES}/${SAMPLE}.2.fq.gz  ${PATH_TO_SAMPLES}/${SAMPLE}_R1_clipped_passed-re-filter.fastq.bz2  ${PATH_TO_SAMPLES}/${SAMPLE}_R2_clipped_passed-re-filter.fastq.bz2
###Keep the .1.fq.gz and .2.fq.gz for the output, as it is needed for stacks to understand that reads are paired when running tsv2bam. 
