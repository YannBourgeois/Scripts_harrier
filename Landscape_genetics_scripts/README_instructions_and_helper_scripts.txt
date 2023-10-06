####The VCF file with which we will work is VCF_M5_N7_80pc.vcf

####Conversion to plink format and preparing input files

vcftools --vcf VCF_M5_N7_80pc.vcf --plink --out VCF_M5_N7_80pc
cut -f1 VCF_M5_N7_80pc.ped > list_in_order_samples.txt

##Below is a little trick to feed plink with more than the standard number of human chromosomes. First column is filled with 1, and SNPs are identified in the 2nd column as contig:position

awk '{print "1", $2, $3, $4}' VCF_M5_N7_80pc.map > VCF_M5_N7_80pc.map2
mv VCF_M5_N7_80pc.map2 VCF_M5_N7_80pc.map
plink --file VCF_M5_N7_80pc --make-bed --out VCF_M5_N7_80pc
plink --file VCF_M5_N7_80pc --recode12 --out VCF_M5_N7_80pc.plink --allow-extra-chr


####Obtain relatedness from NGSRelatev2. Obtain list of unrelated individuals to be used for demographic analyses

./ngsRelate -h VCF_M5_N7_80pc.vcf -O relate.res -T GT -c 1 -z list_in_order_samples.txt -l 0.005
##The -z option outputs the names of the individuals for each pair, which is way more convenient.


####ADMIXTURE (quite fast, but you can also submit a job array to parallelize, cf the stacks folder for examples)

for rep in {1..10}
do
mkdir run_${rep}; cd run_${rep}
cp ../VCF_M5_N7_80pc .
for K in 1 2 3 4 5; 
do 
admixture --cv=10 VCF_M5_N7_80pc.plink.ped $K | tee log${K}.out
done
done

#You can pick the runs with the lowest CV error for each K.
grep -h CV run_*/log*.out

## You can then use a paste command with the list_in_order_samples.txt file to obtain ancestry coefficients with their associated individual. This can be used to produce a file similar to the (provided) ADMIXTURE_K2.txt
paste list_in_order_samples.txt VCF_M5_N7_80pc.plink.2.Q > ADMIXTURE_K2.txt


####For DAPC, we need a specific type of PLINK file

plink --file VCF_M5_N7_80pc --recodeA
mv plink.raw DAPC.raw
cp VCF_M5_N7_80pc.map DAPC.map


####EEMS scripts
##Note: EEMS is rather annoying to install, as it requires outdated libraries AND WILL NOT COMPILE IF USING THE MOST RECENT ONES! Be very careful to follow the instructions at this page. I make it work using conda/mamba (either on my laptop or on a cluster)
##You also need eigen https://gitlab.com/libeigen/eigen/-/releases/3.3.9 
##First create a mamba/conda environment named EEMS. Then you need to install the following libraries:

mamba activate EEMS
gunzip eigen-3.3.9.tar.gz
mamba install -c conda-forge boost=1.74
mamba install -c conda-forge gcc
mamba install -c conda-forge gxx

##You need then to reach the runeems_snps folder and edit the paths in the Makefile (an example below)
EIGEN_INC = /home/yannbourgeois/Desktop/CurrentWORK/Tortoises/EEMS_trial/eigen-3.3.9
BOOST_LIB = /home/yannbourgeois/mambaforge/envs/EEMS/lib/
BOOST_INC = /home/yannbourgeois/mambaforge/envs/EEMS/include

##Then run
make linux
##This should compile, hopefully.

##You can find input files for EEMS (see https://github.com/dipetkov/eems for instructions) in the EEMS folder. From this folder, an EEMS analysis can be launched using the following command:
./runeems_snps --params params-chain1.ini --seed 123
##Here again, don't just do a single run. Create several folders (see the admixture example). Select then the run with the highest final likelihood, or check that they give consistent results.
##The analysis of the ouput is performed in R (see the Quarto document in the folder)


####easySFS + stairwayplot2 scripts


####R scripts to obtain plots and summmary statistics can be found in the Quarto document (R_scripts.qmd)