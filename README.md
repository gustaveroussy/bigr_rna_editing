# rna-editing
This pipeline uses SPRINT and RNAEditingIndexer to identify editing events from paired-end RNA-seq data.

Be careful: this pipeline is in testing!

## Installation
The pipeline is already installed on the Flamingo cluster of Gustave Roussy.  
It is localized here: /mnt/beegfs/pipelines/rna-editing

## Using
You need to make 2 files: a design file and a configuration file.   
### Configuration file
You can copy the example from config/config.yaml.
- **design**: absolute path to your design.csv file.
- **output_dir**: absolute path to the output directory where results will be saved.
- **reference**: the reference to use for the alignment and the idetification of editing events. Possible choices are hg19, hg38, mm10 or mm9. The reference will be downloaded from the UCSC web site.
- **samples_order_for_ggplot**: the order of samples for the x axis of graphs (you can group samples by condition for example). Default is alphabetical order.

### Design file
It must be a comma separated file (.csv where comma is ",") with 3 columns:
- **sample_id**: the sample name of you sample (it could be different that your fastq files).
- **R1_fastq**: absolute path to the R1.fastq.gz file.
- **R2_fastq**: absolute path to the R2.fastq.gz file.

Example:
```
sample_id,R1_fastq,R2_fastq
S1_patient,/mnt/beegfs/scratch/m_aglave/Edition_RNA/data_input/S1-patient_R1.fastq.gz,/mnt/beegfs/scratch/m_aglave/Edition_RNA/data_input/S1-patient_R2.fastq.gz
S2_patient,/mnt/beegfs/scratch/m_aglave/Edition_RNA/data_input/S2-patient_R1.fastq.gz,/mnt/beegfs/scratch/m_aglave/Edition_RNA/data_input/S2-patient_R2.fastq.gz
```
> Notes:
> - sample names mustn't contain special characters or spaces.
> - fastq files must be gzipped.

### Run
You need snakemake (via conda) and singularity (via module load).  
Don't forget to change the path to your configuration file.
```
Example of script:
#!/bin/bash
#using: sbatch run.sh
#SBATCH --job-name=Editing_analysis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=longq

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/m_aglave/.environnement_conda/scRNAseq_10X_user
module load singularity

Editing_pipeline="/mnt/beegfs/pipelines/rna-editing/"
snakemake --profile ${Editing_pipeline}/profiles/slurm \
          -s ${Editing_pipeline}/Snakefile \
          --configfile path_to/my_configuration_file.yaml
```

## Steps of the pipeline
1. Symbolic link of fastq files
2. QC & Trimming (fastQC, fastp & multiqc)
3. BWA index generation (via SPRINT)
3. BWA alignement (via SPRINT)
4. Identification of Editing events (SPRINT)
5. Summary of SPRINT results (R)
6. Bam sorting (Samtools)
7. Identification of Editing events (RNAEditingIndexer)
8. Summary of RNAEditingIndexer results (R)

Information about Editing tools:  
SPRINT:  
https://github.com/jumphone/SPRINT  
https://academic.oup.com/bioinformatics/article/33/22/3538/4004872  
RNAEditingIndexer:  
https://github.com/a2iEditing/RNAEditingIndexer  
https://pubmed.ncbi.nlm.nih.gov/31636457/  

