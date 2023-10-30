#!/bin/bash
#SBATCH --mem=32G
#SBATCH -p sciama3.q  ####change the queue name to match your own cluster.
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --array=1-20        # 5 individuals times 4 intervals
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


config=list_fastq_jobs_rnaseq_gatk.txt

SAMPLE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
CHUNK=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

gatk HaplotypeCaller \
-R GCF_929443795.1_bAccGen1.1_genomic.fna \
-I RGsplitN_mrkdup_${SAMPLE}.bam \
-ERC GVCF \
-O ${SAMPLE}_interval_${CHUNK}.g.vcf.gz \
-dont-use-soft-clipped-bases \
-L intervals${CHUNK}.list \
--standard-min-confidence-threshold-for-calling 20 


