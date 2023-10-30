######We can download individuals from the close relative C. melanoleucos using the prefetch command from sra-tools.
######See https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump for a tutorial.
######First we download the sister species data  (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3203217)

prefetch SRR3203217
fasterq-dump --split-files SRR3203217/SRR3203217.sra


###Prepare a file for the arrays as follows:

ls *_1.fastq.gz -d | sed "s:_1.fastq.gz::" | nl -w2  >> list_fastq_files_rnaseq.txt

###This should give a file looking like this:

1	CDS2347
2	SRR3203217
3	D252517
4	DA252518
5	DA252521
6	DA282337


#####Align reads on the reference genome using STAR

##First, we trim the reads (see associated scripts for more details):
sbatch trimming_array.sh

##We prepare the genome index for STAR
sbatch prepare_genome_ref_for_star.sh

##We run the alignment + prepare the resulting BAM files for genotyping with GATK. Script runs STAR, MarkDuplicates, SplitNCigarReads and AddOrReplaceReadGroups
sbatch star_array_mark_duplicates.sh

#####Variant calling using GATK (following protocol suggested at https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
###We split the genome into four chunks (coordinates in files named intervals${CHUNK}.list) and parallelize the calling across the 5 individuals for the four chunks. This makes 20 jobs to run.
###We prepare the config file for the array, for example:

ls *_1.fastq.gz -d | sed "s:_1.fastq.gz:\t 1:"  >> list_fastq_jobs_rnaseq_gatk.txt
ls *_1.fastq.gz -d | sed "s:_1.fastq.gz:\t 2:"  >> list_fastq_jobs_rnaseq_gatk.txt
ls *_1.fastq.gz -d | sed "s:_1.fastq.gz:\t 3:"  >> list_fastq_jobs_rnaseq_gatk.txt
ls *_1.fastq.gz -d | sed "s:_1.fastq.gz:\t 4:"  >> list_fastq_jobs_rnaseq_gatk.txt

cat list_fastq_jobs_rnaseq_gatk.txt | nl -w2 > list_fastq_jobs_rnaseq_gatk.txt2;mv list_fastq_jobs_rnaseq_gatk.txt2 list_fastq_jobs_rnaseq_gatk.txt

##Launching the calling.
sbatch gatk_haplocaller_array.sh



###Joint genotyping
###See https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport
###this option is important to speed up the calling with many contigs, we keep it even with gentilis genome which is less fragmented than the initial accipiter nisus.
### --merge-contigs-into-num-partitions
###See also https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines

sbatch joint_genotyping_array.sh



######We can then merge the VCF files for the four intervals as follows:

gatk MergeVcfs \
          I=fivesamples_interval_1.vcf.gz \
          I=fivesamples_interval_2.vcf.gz \
          I=fivesamples_interval_3.vcf.gz \
          I=fivesamples_interval_4.vcf.gz \
          O=all_samples.vcf.gz


######We can extract summary statistics for quality with an interactive run:

srun --mem 3GB --time 5:00:00 --pty /bin/bash

bcftools query fivesamples_interval_1.vcf.gz -f '%QUAL\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > qual_stats.QUAL.FS.SOR.MQRS.RPRS.QD.MQ.DP.txt

###See http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/

#FS (GATK recommendation FS > 60): FS > 60 (seems fine given the distribution)
#SOR (GATK recommendation SOR > 3): SOR > 3 (this seems just fine)
#MQ (GATK recommendation MQ < 40): MQ < 40 (we would loose many variants with a lower threshold, so we stick to the recommended)
#MQRankSum (GATK recommendation MQRankSum > -12.5): MQRankSum > -5.0 and MQRankSum < 5.0 (a value around 0 is preferred and we don’t loose many variants if we set strict thresholds)
#QD (GATK recommendation QD < 2): QD < 2 fine. Very high quality data. We can afford it.
#ReadPosRankSum (GATK recommendation ReadPosRankSum > -8.0): ReadPosRankSum < -4.0 (not many variants are below -4, so we can set it more strict)
#INFO/DP (no GATK recommendation, dependend on sample number and sequencing depth): INFO/DP > 500 (the cumulated DP is roughly 100X, we exclude sites that show about twice or half the mean DP)

###We also exclude genotypes that have less than GQ 20, depth of 8X.

bcftools filter -e 'FS>60.0 || SOR>3 || MQ<40 || MQRankSum<-5.0 || MQRankSum>5.0 || QD<2.0 || ReadPosRankSum<-4.0 || INFO/DP>101' all_samples.vcf.gz -O z -o all_samples_filtered_site_level.vcf.gz 
  
vcftools --gzvcf all_samples_filtered_site_level.vcf.gz --max-missing-count 2 --mac 1 --minGQ 20 --minDP 8 --max-mean_DP 100 --out filtered_vcf_all_samples --recode --recode-INFO-all
###Discovered that the --max-missing-count options counts the NUMBER OF ALLELES!!!! in VCFTOOLS 1.16. Here we have at most one missing individual


####Annotation. We can use SNPEff. Originally we used Accipiter nisus which was prebuilt. Accipiter gentilis is not yet in there, so we need to follow the instructions at https://pcingola.github.io/SnpEff/snpeff/build_db_gff_gtf/#step-2-build-using-gene-annotations-and-reference-sequences

snpEff Accipiter_nisus_ver1.0.99 filtered_vcf_all_samples.recode.vcf -V > all.ann.vcf




