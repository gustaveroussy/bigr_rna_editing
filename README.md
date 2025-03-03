# bigr_rna_editing
This pipeline uses SPRINT and RNAEditingIndexer to identify editing events from paired-end RNA-seq data.

## Installation

### Gustave Roussy users
The pipeline is already installed on the Flamingo cluster of Gustave Roussy.  
It is localized here: /mnt/beegfs/pipelines/bigr_rna_editing/<version>

### Admin : installation of a new version of the pipeline
#### :one: Download pipeline
```
cd /mnt/beegfs/pipelines/bigr_rna_editing/
VERSION="1.1.0"
git clone https://github.com/gustaveroussy/bigr_rna_editing.git ${VERSION}
```
#### :two: Download Singularity images
Download all singularity images from [Zenodo](https://zenodo.org/records/14916660):
```
cd /mnt/beegfs02/pipelines/bigr_rna_editing/${VERSION}/envs/singularity/
wget https://zenodo.org/api/records/14916660/files-archive
unzip files-archive
```
#### :three: Install Snakemake environment
```
source /mnt/beegfs02/software/recherche/miniconda/25.1.1/etc/profile.d/conda.sh
conda env create -f /mnt/beegfs02/pipelines/bigr_rna_editing/${VERSION}/envs/conda/snakemake.yaml --prefix=/mnt/beegfs02/pipelines/bigr_rna_editing/${VERSION}/envs/compiled_conda/snakemake -y
```
You are now ready to use the pipeline!

## Using
You need to make 2 files: a design file and a configuration file.   
### Configuration file
- **design**: absolute path to your design.csv file.
- **output_dir**: absolute path to the output directory where results will be saved.
- **reference**: the reference to use for the alignment and the idetification of editing events. Possible choices are "hg19", "hg38", "mm10" or "mm9". The reference will be downloaded from the UCSC web site.
- **samples_order_for_ggplot** (optional): the order of samples for the x axis of graphs (you can order samples by condition for example). Default is alphabetical order.
- **SPRINT_extra** (optional): extra parameters for "SPRINT main" command.
- **RNAEditingIndexer_extra** (optional): extra parameters for "RNAEditingIndexer" command.
- **nb_sampled_reads** (optional): number of reads to sample. Possible choices are "" for no sampling, "50000000" for 50M of reads (can be another integer), or "auto" (to sample the minimum number of reads obtain throught all samples if they have more than _min_nb_sampled_reads_for_auto_ option, else the threshold is the value of _min_nb_sampled_reads_for_auto_. If a sample has less than the threshold, all its reads are used.
- **min_nb_sampled_reads_for_auto** (optional): minimum number of reads to sample if _nb_sampled_reads_ is set to "auto". Defaut is 50M of reads.

Example:
```
design: "/mnt/beegfs02/scratch/m_aglave/Editing_analysis/script/design.csv"
output_dir: "/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_output/"
reference: "hg38"
samples_order_for_ggplot: "S1_patient,S3_patient,S2_patient"
SPRINT_extra: ""
RNAEditingIndexer_extra: ""
```
### Design file
It must be a comma separated file (.csv where comma is ",") with 3 columns:
- **sample_id**: the sample name of you sample (it could be different that your fastq files).
- **R1_fastq**: absolute path to the R1.fastq.gz file.
- **R2_fastq**: absolute path to the R2.fastq.gz file.

Example:
```
sample_id,R1_fastq,R2_fastq
S1_patient,/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_input/S1-patient_R1.fastq.gz,/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_input/S1-patient_R2.fastq.gz
S2_patient,/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_input/S2-patient_R1.fastq.gz,/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_input/S2-patient_R2.fastq.gz
S3_patient,/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_input/S3-patient_R1.fastq.gz,/mnt/beegfs02/scratch/m_aglave/Editing_analysis/data_input/S3-patient_R2.fastq.gz
```
> Notes:
> - sample names mustn't contain special characters or spaces.
> - fastq files must be gzipped.

### Run
You need snakemake (via conda) and singularity (via module load). They are already installed for you on Flamingo, just follow the example below.  
Don't forget to change the version of the pipeline and the path to your configuration file.  
Example of script:
```
#!/bin/bash
#using: sbatch run.sh
#SBATCH --job-name=Editing_analysis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250M
#SBATCH --partition=longq

source /mnt/beegfs02/software/recherche/miniconda/25.1.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs02/pipelines/bigr_rna_editing/<version>/envs/conda/snakemake
module load singularity
Editing_pipeline="/mnt/beegfs02/pipelines/bigr_rna_editing/<version>/"

snakemake --profile ${Editing_pipeline}/profiles/slurm \
          -s ${Editing_pipeline}/Snakefile \
          --configfile <path_to/my_configuration_file.yaml>
```

## Steps of the pipeline
1. Symbolic link of fastq files
2. Reads QC & Trimming (fastQC, fastp & multiqc)
3. BWA index generation (via SPRINT)
4. Sampling read (optional) (seqtk)
5. BWA alignement (via SPRINT)
6. Aligment QC (Samtools & multiqc)
7. Identification of Editing events (SPRINT) (this step takes 2-3 days!)
8. Summary of SPRINT results (R)
9. Bam sorting (Samtools)
10. Identification of Editing events (RNAEditingIndexer)
11. Summary of RNAEditingIndexer results (R)

Information about Editing tools:  
SPRINT:  
https://github.com/jumphone/SPRINT  
https://academic.oup.com/bioinformatics/article/33/22/3538/4004872  
RNAEditingIndexer:  
https://github.com/a2iEditing/RNAEditingIndexer  
https://pubmed.ncbi.nlm.nih.gov/31636457/  
