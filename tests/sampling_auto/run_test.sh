#!/bin/bash

########################################################################
## Script to launch RNAseq Editing pipeline
##
## using: sbatch /mnt/beegfs/pipelines/bigr_rna_editing/dev/tests/sampling_auto/run_test.sh
##
########################################################################

## JOB PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#SBATCH --job-name=Editing_tests_sampling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=bigmemq

source /mnt/beegfs/software/miniconda/24.3.0/etc/profile.d/conda.sh
conda activate /mnt/beegfs/pipelines/bigr_rna_editing/dev/envs/conda/snakemake
module load singularity
Editing_pipeline="/mnt/beegfs/pipelines/bigr_rna_editing/dev/"

snakemake --profile ${Editing_pipeline}/profiles/slurm \
          -s ${Editing_pipeline}/Snakefile \
          --configfile /mnt/beegfs/pipelines/bigr_rna_editing/dev/tests/sampling_auto/config.yaml
          #--attempt 2 --restart-times 5