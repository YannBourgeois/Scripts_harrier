#### You will find in this repository job arrays (for SLURM schedulers) to run the stacks pipeline for multiple values of N and M
####You should run the scripts in this order: cutadapt.sh > ustacks_array.sh > stacks_last_steps.sh
####This README file lists the intermediate commands that you need to run to make sure the scripts will run properly. Read it carefully!
#### For more information about slack and its options, please consult https://catchenlab.life.illinois.edu/stacks/manual/
#### It should be possible to create a snakemake pipeline to run all these scripts in an automated way, but at the moment you need to submit the job arrays corresponding to the different steps yourself.

#### We first need to generate the configuration files that will list our samples, the value of mismatches M that we want to use to group reads into a contig for each individual
#### The list_fastq_files.txt file lists the array number, the ID number (for stacks), the Name of the sample, and the M number (can be 2, 5, 7 and 9)
#### In our case, the raw fastq files are named as follows: /mnt/lustre/yb24/Harrier_analyses/RE_processed/P01-C02-11_R1_clipped_passed-re-filter.fastq.bz2
#### You can change this by whatever suffix works with your own data of course.
#### Below is a possible, not super elegant, semi-automated way to generate the file that will be read when launching an array job.
#### Replace the path in PATH_TO_SAMPLES by your own path to the analysis folder

PATH_TO_SAMPLES=/mnt/lustre/yb24/Harrier_analyses/RE_processed/
for m in 1 3 5 7
do
ls *_R1_clipped_passed-re-filter.fastq.bz2 -d | sed "s:\$PATH_TO_SAMPLES::" | sed "s:_R1_clipped_passed-re-filter.fastq.bz2:\t${m}:" | nl -w2  >> list_fastq_files.txt
done
cat list_fastq_files.txt | nl -w2 > list_fastq_files.txt2;mv list_fastq_files.txt2 list_fastq_files.txt

#### We only need to run cutadapt once on each sample

head -n61 list_fastq_files.txt > list_for_cutadapt.txt

####Create the folders in which the ustacks output will be placed
mkdir ${PATH_TO_SAMPLES}stacks_M1
mkdir ${PATH_TO_SAMPLES}stacks_M3
mkdir ${PATH_TO_SAMPLES}stacks_M5
mkdir ${PATH_TO_SAMPLES}stacks_M7


##### You can then run the script cutadapt.sh, then ustacks_array.sh, using sbatch followed by the script's name. Have a look at the scripts to understand what it does and edit the paths to the work directories to match your own folder organization
##### Once these two scripts have run, we generate a population file listing individuals
##### It can be easily obtained as follows, assuming you are in the initial working directory (here /mnt/lustre/yb24/Harrier_analyses/RE_processed/):
cd stacks_M1
ls -s *tags*  | sed "s:.tags.tsv.gz:\tPop1:" > ../populations_cstacks.txt
cd ..


#####We also generate output folders for different N (cstacks) and M (ustacks). Generate symbolic links to the initial ustacks folders.
for m in 1 3 5 7
do
n=${m}
mkdir stacks_M${m}_N${n}
cd stacks_M${m}_N${n}
ln -s ../stacks_M${m}/*1.alleles.tsv.gz .
ln -s ../stacks_M${m}/*1.snps.tsv.gz .
ln -s ../stacks_M${m}/*1.tags.tsv.gz .
cd ..

n=`expr $m + 2`
mkdir stacks_M${m}_N${n}
cd stacks_M${m}_N${n}
ln -s ../stacks_M${m}/*1.alleles.tsv.gz .
ln -s ../stacks_M${m}/*1.snps.tsv.gz .
ln -s ../stacks_M${m}/*1.tags.tsv.gz .
cd ..

if [ $m -ge 3 ]
then
n=`expr $m - 2`
mkdir stacks_M${m}_N${n}
cd stacks_M${m}_N${n}
ln -s ../stacks_M${m}/*1.alleles.tsv.gz .
ln -s ../stacks_M${m}/*1.snps.tsv.gz .
ln -s ../stacks_M${m}/*1.tags.tsv.gz .
cd ..

fi
done

####At last, we generate a list of combinations for the job array to run the last steps of the stacks pipeline 

ls -d stacks_M*_N* | nl -w2 > list_combinations_folders_cstacks.txt
ls -d stacks_M*_N* | cut -f3 -d "_" | sed "s:N::" | paste list_combinations_folders_cstacks.txt - > list_combinations_folders_cstacks.txt2;mv list_combinations_folders_cstacks.txt2 list_combinations_folders_cstacks.txt

####You can now run the stacks_last_steps.sh script
####Once all is done, you can use these simple commands using vcftools to quickly check the set of parameters producing the most loci and SNPs
####Mean + sd coverage of assembled loci

grep cov stacks_M*_N*/*log

####Number of polymorphic loci and SNPs

for i in stacks_M*_N*
do
vcftools --vcf $i/populations.snps.vcf --mac 1 --max-missing 0.8 --minGQ 20 --minDP 6 --recode --hwe 0.0001 --recode-INFO-all --remove-indv P02-D03-59_R1 --remove-indels --max-alleles 2 --out stacks_filtered_80pc
echo $i   ###Which combination
cut -f1 stacks_filtered_80pc.recode.vcf | grep -v "#" | sort | uniq | wc -l  ###How many polymorphic contigs
grep -v "#" stacks_filtered_80pc.recode.vcf -c    ###How many variants
done



