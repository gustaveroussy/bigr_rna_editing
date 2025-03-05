#!/bin/bash

########################################################################
## Script to launch RNAseq Editing pipeline
##
## using: sbatch /mnt/beegfs02/pipelines/bigr_rna_editing/1.1.1/tests/sampling_auto/run_test.sh
##
########################################################################

## JOB PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#SBATCH --job-name=Editing_tests_sampling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=longq

source /mnt/beegfs02/software/recherche/miniconda/25.1.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs02/pipelines/bigr_rna_editing/1.1.1/envs/compiled_conda/snakemake
module load singularity-ce
Editing_pipeline="/mnt/beegfs02/pipelines/bigr_rna_editing/1.1.1/"

snakemake --profile ${Editing_pipeline}/profiles/slurm \
          -s ${Editing_pipeline}/Snakefile \
          --configfile ${Editing_pipeline}/tests/sampling_auto/config.yaml
