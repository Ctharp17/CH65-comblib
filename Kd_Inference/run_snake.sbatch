#!/bin/bash
#
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J mergedKdv2 # Job name
#SBATCH -n 10 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 100000 # Memory in MB
#SBATCH -o /n/desai_lab/users/aphillips/babyYODA/titeseq/CH65_full_titeseq_analysis_v7/Kd_inference/slurm_outputs/snakemake_%A_%a.out
#SBATCH -e /n/desai_lab/users/aphillips/babyYODA/titeseq/CH65_full_titeseq_analysis_v7/Kd_inference/slurm_outputs/snakemake_%A_%a.err

module load python/3.7.7-fasrc01
source activate omicron

snakemake -j 10
scp -r results_CH65 ../../../../../../../desai_lab/users/aphillips/babyYODA/titeseq/CH65_full_titeseq_analysis_v7/Kd_inference/
