#!/bin/bash
#SBATCH --mem=15G
#SBATCH -p sciama3.q
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --array=1-5        # 5 individuals 
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Walltime format hh:mm:ss
#SBATCH --time=24:00:00
# Output and error files
#SBATCH -o job.%J.out
#SBATCH -e job.%J.err
# **** Put all #SBATCH directives above this line! ****
# **** Otherwise they will not be in effective! ****
#
# **** Actual commands start here ****
# Load modules here (safety measure)
config=list_fastq_files_rnaseq.txt

SAMPLE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


trimmomatic PE ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz trimmed_ ${SAMPLE}_F.fq.gz 
trimmed_unpaired${SAMPLE}_F.fq.gz trimmed_${SAMPLE}_R.fq.gz 
trimmed_unpaired${SAMPLE}_R.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:TRUE LEADING:3 TRAILING:3 MINLEN:30 SLIDINGWINDOW:4:20

###Generate FASTQC reports
fastqc trimmed_${SAMPLE}_F.fq.gz
fastqc trimmed_${SAMPLE}_R.fq.gz

