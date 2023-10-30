#!/bin/bash
#SBATCH --mem=32G
#SBATCH -p sciama3.q
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --array=1-5        # 5 individuals 
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=8
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

STAR \
--genomeDir genome_accipiter \
--runThreadN 8 \
--readFilesIn trimmed_${SAMPLE}_F.fq.gz trimmed_${SAMPLE}_R.fq.gz \
--readFilesCommand "gunzip -c" \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--outFileNamePrefix ${SAMPLE} --outFilterMismatchNmax 15 


gatk MarkDuplicates \
--INPUT ${SAMPLE}Aligned.sortedByCoord.out.bam \
--OUTPUT mrkdup_${SAMPLE}.bam  \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY SILENT \
--METRICS_FILE ${SAMPLE}.metrics

gatk SplitNCigarReads \
-R GCF_929443795.1_bAccGen1.1_genomic.fna \
-I mrkdup_${SAMPLE}.bam \
-O splitN_mrkdup_${SAMPLE}.bam 

gatk AddOrReplaceReadGroups \
I=splitN_mrkdup_${SAMPLE}.bam \
O=RGsplitN_mrkdup_${SAMPLE}.bam \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=${SAMPLE} CREATE_INDEX=true

