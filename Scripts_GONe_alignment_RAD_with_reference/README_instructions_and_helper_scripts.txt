#### First we index the reference genome
####Accipiter gentilis is a more complete assembly, which is what matters most.
####https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/929/443/795/GCA_929443795.1_bAccGen1.1/GCA_929443795.1_bAccGen1.1_genomic.fna.gz

bwa index GCA_929443795.1_bAccGen1.1_genomic.fna

###Aligning the reads on the reference:
sbatch bwa_commands.sh

####List all resulting bam files in a file:
ls *bam > file_bam.list

###SNP calling and filtering with freebayes/vcftools:
sbatch freebayes.sh


####This is an annoying step, but has to be done to mage GONe work (requires numbers as chromosome IDs).

sed -i -e "s:0\tOV839361.1:1\tOV839361.1:" trial*.map
sed -i -e "s:0\tOV839370.1:2\tOV839370.1:" trial*.map
sed -i -e "s:0\tOV839372.1:3\tOV839372.1:" trial*.map
sed -i -e "s:0\tOV839373.1:4\tOV839373.1:" trial*.map
sed -i -e "s:0\tOV839374.1:5\tOV839374.1:" trial*.map
sed -i -e "s:0\tOV839375.1:6\tOV839375.1:" trial*.map
sed -i -e "s:0\tOV839376.1:7\tOV839376.1:" trial*.map
sed -i -e "s:0\tOV839377.1:8\tOV839377.1:" trial*.map
sed -i -e "s:0\tOV839378.1:9\tOV839378.1:" trial*.map
sed -i -e "s:0\tOV839379.1:10\tOV839379.1:" trial*.map
sed -i -e "s:0\tOV839380.1:11\tOV839380.1:" trial*.map
sed -i -e "s:0\tOV839362.1:12\tOV839362.1:" trial*.map
sed -i -e "s:0\tOV839381.1:13\tOV839381.1:" trial*.map
sed -i -e "s:0\tOV839382.1:14\tOV839382.1:" trial*.map
sed -i -e "s:0\tOV839383.1:15\tOV839383.1:" trial*.map
sed -i -e "s:0\tOV839384.1:16\tOV839384.1:" trial*.map
sed -i -e "s:0\tOV839385.1:17\tOV839385.1:" trial*.map
sed -i -e "s:0\tOV839386.1:18\tOV839386.1:" trial*.map
sed -i -e "s:0\tOV839387.1:19\tOV839387.1:" trial*.map
sed -i -e "s:0\tOV839388.1:20\tOV839388.1:" trial*.map
sed -i -e "s:0\tOV839389.1:21\tOV839389.1:" trial*.map
sed -i -e "s:0\tOV839390.1:22\tOV839390.1:" trial*.map
sed -i -e "s:0\tOV839363.1:23\tOV839363.1:" trial*.map
sed -i -e "s:0\tOV839391.1:24\tOV839391.1:" trial*.map
sed -i -e "s:0\tOV839392.1:25\tOV839392.1:" trial*.map
sed -i -e "s:0\tOV839393.1:26\tOV839393.1:" trial*.map
sed -i -e "s:0\tOV839394.1:27\tOV839394.1:" trial*.map
sed -i -e "s:0\tOV839395.1:28\tOV839395.1:" trial*.map
sed -i -e "s:0\tOV839396.1:29\tOV839396.1:" trial*.map
sed -i -e "s:0\tOV839397.1:30\tOV839397.1:" trial*.map
sed -i -e "s:0\tOV839398.1:31\tOV839398.1:" trial*.map
sed -i -e "s:0\tOV839399.1:32\tOV839399.1:" trial*.map
sed -i -e "s:0\tOV839364.1:33\tOV839364.1:" trial*.map
sed -i -e "s:0\tOV839365.1:34\tOV839365.1:" trial*.map
sed -i -e "s:0\tOV839366.1:35\tOV839366.1:" trial*.map
sed -i -e "s:0\tOV839367.1:36\tOV839367.1:" trial*.map
sed -i -e "s:0\tOV839368.1:37\tOV839368.1:" trial*.map
sed -i -e "s:0\tOV839369.1:38\tOV839369.1:" trial*.map
sed -i -e "s:\t: :g" trial*.*
sed -i -e 's/R1 0 0 0 0/R1 0 0 0 -9 /' trial*.ped  ###No sex indication. Sex chromosome is excluded from analyses 


###GONe extracts the number of chromosomes by looking at the number starting the last row, so that needs to be carefully ordered
###The scripts can be obtained at https://github.com/esrud/GONE/tree/master/Linux
####2cM/Mb seems like a reasonable value. https://avianres.biomedcentral.com/articles/10.1186/s40657-018-0096-7/tables/1
####We edit the INPUT_PARAMETERS_FILE accordingly.

./script_GONE.sh trial_popall
./script_GONE.sh trial_North
./script_GONE.sh trial_South
