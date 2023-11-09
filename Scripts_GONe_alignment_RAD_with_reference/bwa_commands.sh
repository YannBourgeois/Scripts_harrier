#!/bin/bash
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --job-name=bwa_job
#SBATCH --array=1-61         # 61 individuals
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=12  ###Better if this matches the -t parameter in bwa
#SBATCH --output=bwa_job.%J.out
#SBATCH --error=bwa_job.%J.err
## ***Commands starts here***

config=list_for_cutadapt.txt ###We can reuse the same fastq as the ones from stacks. Third column is the sample.
PATH_TO_SAMPLES=/mnt/lustre/yb24/Harrier_analyses/RE_processed/Second_reads_trimmed/  ###Edit so it points to the folder with the fastq files obtained with cutadapt
SAMPLE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config) 
bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" -M -t 12 GCA_929443795.1_bAccGen1.1_genomic.fna ${PATH_TO_SAMPLES}/${SAMPLE}.1.fq.gz ${PATH_TO_SAMPLES}/${SAMPLE}.2.fq.gz | samtools sort -@12 -o ${SAMPLE}.bam -
