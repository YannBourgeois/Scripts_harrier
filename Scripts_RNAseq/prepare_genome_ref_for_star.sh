#!/bin/bash
#SBATCH --mem=32G
#SBATCH -p sciama3.q ###local queue name. change this depending on your cluster
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
##
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
STAR \
--runMode genomeGenerate \
--genomeDir genome_accipiter \
--genomeFastaFiles GCF_929443795.1_bAccGen1.1_genomic.fna \
--sjdbGTFfile GCF_929443795.1_bAccGen1.1_genomic.gtf \
--sjdbOverhang 125 \
--runThreadN 8 \
--limitGenomeGenerateRAM=32000000000 --genomeChrBinNbits=14