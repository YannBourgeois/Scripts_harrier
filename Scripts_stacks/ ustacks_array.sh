#!/bin/bash
#SBATCH --mem=12GB
#SBATCH --nodes=1
#SBATCH --job-name=stacks_job
#SBATCH --array=1-244         # 61 individuals times 4 values for M jobs
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4  ###Better if this matches the -p parameter in ustacks
#SBATCH --output=stacks_job.%J.out
#SBATCH --error=stacks_job.%J.err
## ***Commands starts here***

####This script will run ustacks with M values of 1, 3, 5 and 7 (check the config file) for each of the 61 samples, and place the output in folders named stacks_M1, stacks_M3, etc.
config=list_fastq_files.txt

PATH_TO_SAMPLES=/mnt/lustre/yb24/Harrier_analyses/RE_processed/Second_reads_trimmed/
ID_NUMBER=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
SAMPLE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
M_VALUE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

ustacks -t gzfastq -f ${PATH_TO_SAMPLES}/${SAMPLE}.1.fq.gz -o ./stacks_M${M_VALUE} -i ${ID_NUMBER} -m 3 -M ${M_VALUE} -p 4