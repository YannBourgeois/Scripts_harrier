!/bin/bash
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --job-name=stacks_job
#SBATCH --array=1-11        # 11 combinations of M and N
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=10   #Might be too much, if lowering this number, make sure to also edit the -p or -t options below
#SBATCH --output=stacks_job.%J.out
#SBATCH --error=stacks_job.%J.err
## ***Commands starts here***

####This script will run ustacks with M values of 1, 3, 5 and 7 (check the config file) for each of the 61 samples, and place the output in folders named stacks_M1, stacks_M3, etc.
config=list_combinations_folders_cstacks.txt

PATH_TO_SAMPLES=/mnt/lustre/yb24/Harrier_analyses/RE_processed/Second_reads_trimmed/
FOLDER=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
N_VALUE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)


cstacks -P ${FOLDER} -M ./populations_cstacks.txt -n ${N_VALUE} -p 10
sstacks -P ${FOLDER} -M populations_cstacks.txt -p 10
tsv2bam -P ${FOLDER} -R ${PATH_TO_SAMPLES} -M populations_cstacks.txt -t 10
gstacks -P ${FOLDER} -M populations_cstacks.txt -t 10
