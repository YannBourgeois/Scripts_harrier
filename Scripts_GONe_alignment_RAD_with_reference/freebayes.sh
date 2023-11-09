#!/bin/bash
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --job-name=freebayes_job
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1  ###Better if this matches the -t parameter in bwa
#SBATCH --output=freebayes_job.%J.out
#SBATCH --error=freebayes_job.%J.err
## ***Commands starts here***

freebayes -L file_bam.list -f GCA_929443795.1_bAccGen1.1_genomic.fna --genotype-qualities | vcfallelicprimitives -kg > freebayes_RADseq_raw
vcftools --vcf freebayes_RADseq_raw.vcf --mac 1 --minQ 20 --minGQ 20 --minDP 6 --max-alleles 2 --max-missing 0.8 --remove-indv P02-D03-59_R1 --remove-indels --hardy --out HWExcess_filter

###removing putative paralogs with heterozygotes excess. Decided on this threshold based on empirirical distributions, should be changed depending on what p-value corresponds to clear issue (i.e. almost all individuals are heterozygous).

gawk '{ if ($8 <= 0.0001) { print } }' HWExcess_filter.hwe | cut -f1,2 > HWExcess_filter.txt 

vcftools --vcf freebayes_RADseq_raw.vcf --mac 1 --minQ 20 --minGQ 20 --minDP 6 --max-alleles 2 --exclude-positions HWExcess_filter.txt --max-missing 0.8 --remove-indv P02-D03-59_R1 --remove-indels --plink --out trial_popall
vcftools --vcf freebayes_RADseq_raw.vcf --mac 1 --minQ 20 --minGQ 20 --minDP 6 --max-alleles 2 --exclude-positions HWExcess_filter.txt --max-missing 0.8 --remove-indv P02-D03-59_R1 --remove-indels --recode --recode-INFO-all --out trial_popall
vcftools --vcf trial_popall.recode.vcf --mac 1 --keep North --max-missing 0.8 --plink --out trial_North
vcftools --vcf trial_popall.recode.vcf --mac 1 --keep South --max-missing 0.8 --plink --out trial_South
####North and South contain lists of individuals belonging to the North-Eastern and South-Western groups.
