#!/bin/bash
#SBATCH --mail-user=lawrence.fung@sjsu.edu --mail-type=END,FAIL
#SBATCH --job-name=nf_ref_bowtie_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000

# module load conda/conda
# # activate the conda environment (use source and not just $conda activate)
# try without # SBATCH --time=10:00:00

# source /scratch/home/lfung/miniconda3/envs/bin/activate

start=$(date +%s)
echo "start time: $start" # Print the timestamp in seconds

echo "hostname: $HOSTNAME"


# automatically detects where Conda is installed without needing a hardcoded source path
eval "$(conda shell.bash hook)"
echo "Activating Conda Environment"
conda activate nf_code

###########################
###########################

cd '/scratch/home/lfung/PFAS_Data_NF_ref_bowtie'

# srun hostname # print the name of the node
nextflow run funghub/PFAS_project --input PRJNA1137368 -profile spartan_hpc -latest

###########################
###########################

conda deactivate

end=$(date +%s)

echo "end time: $end" # Print the timestamp in seconds

runtime=$((end - start))

echo "run time: $runtime" # print the run time

# run UNIX command hostname on the node, print it in output
srun hostname

# wait for 60 secs then exit job step
srun sleep 60

# sbatch this file to run nextflow