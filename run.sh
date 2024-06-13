#!/bin/bash

#using: sbatch /mnt/beegfs/scratch/m_aglave/test_Edition_TCGA/script/run.sh

#SBATCH --job-name=test_Editing_pipeline
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=longq

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/m_aglave/.environnement_conda/scRNAseq_10X_user
module load singularity

snakemake --profile /mnt/beegfs/scratch/m_aglave/test_Edition_TCGA/script/profiles/slurm -s /mnt/beegfs/scratch/m_aglave/test_Edition_TCGA/script/Snakefile --configfile /mnt/beegfs/scratch/m_aglave/test_Edition_TCGA/script/config/config.yaml
