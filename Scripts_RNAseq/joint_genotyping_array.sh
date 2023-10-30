#!/bin/bash
#SBATCH --mem=32G
#SBATCH -p sciama3.q  ####change the queue name to match your own cluster.
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --array=1-4       # 4 intervals
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



gatk GenomicsDBImport \
-V CDS2347_interval_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
-V SRR3203217_interval_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
-V D252517_interval_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
-V DA252518_interval_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
-V DA252521_interval_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
-V DA282337_interval_${SLURM_ARRAY_TASK_ID}.g.vcf.gz \
  --genomicsdb-workspace-path ./interval_${SLURM_ARRAY_TASK_ID}_db \
  --tmp-dir /mnt/lustre/yb24/tmp \    ####edit so it matches your own workspace
--merge-input-intervals true --merge-contigs-into-num-partitions 30 -L intervals${SLURM_ARRAY_TASK_ID}.list


TILEDB_DISABLE_FILE_LOCKING=1 gatk GenotypeGVCFs \   ####Remove the TILEDB_DISABLE_FILE_LOCKING=1 if your cluster does not use a LUSTRE file system
   -R GCF_929443795.1_bAccGen1.1_genomic.fna \
   -V gendb://interval_${SLURM_ARRAY_TASK_ID}_db \
   -O fivesamples_interval_${SLURM_ARRAY_TASK_ID}.vcf.gz