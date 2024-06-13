import pandas as pd
import os

sys.stderr.write("\n############################################################ \n\n")

sys.stderr.write("                           ++++++               \n")
sys.stderr.write("                        ++++++++++++            \n")
sys.stderr.write("                     ++++++++++++++++++         \n")
sys.stderr.write("                  +++++++++++++++++++++++++     \n")
sys.stderr.write("               ++++++++++++++++++++++++++++++   \n")
sys.stderr.write("               ++++|  _ \(_)/ ____|  __ \++++   \n")
sys.stderr.write("               ++++| |+) |_| |++++| |++) |+++   \n")
sys.stderr.write("               ++++|  _ <| | |+|_ |  _  /++++   \n")
sys.stderr.write("               ++++| |+) | | |++| | |+\ \++++   \n")
sys.stderr.write("               ++++|____/|_|\_____|_|++\_\+++   \n")
sys.stderr.write("               ++++++++++++++++++++++++++++++   \n")
sys.stderr.write("                  +++++++++++++++++++++++++     \n")
sys.stderr.write("                     +++++++++++++++++++        \n")
sys.stderr.write("                        +++++++++++++           \n")
sys.stderr.write("                           +++++++              \n")

sys.stderr.write("\n                   Editing RNA pipeline       \n")

sys.stderr.write("\nFor any question, sent an email to bigr@gustaveroussy.fr")
sys.stderr.write("\n############################################################ \n\n")

### parameters ###################################################################################################################################

#pipeline directory
PIPELINE_DIR = workflow.snakefile
PIPELINE_DIR = PIPELINE_DIR.replace("/Snakefile", "")

sys.stderr.write("\n################### Checking Config and Design File ###################\n\n")
#output directory
if "output_dir" in config:
    OUTPUT_DIR = config["output_dir"]
    sys.stderr.write("Output directory set to " + OUTPUT_DIR + ".\n")
else:
    OUTPUT_DIR = os.getcwd()
    sys.stderr.write("No 'output_dir' found in config file, output directory set to " + OUTPUT_DIR + ".\n")

#design file
if "design" in config :
    if os.path.isfile(config["design"]):
        sys.stderr.write("Design file found: " + config["design"] + ".\n")
    else: 
        sys.exit("Error: 'design' file not found at the given path. Check your configfile.\n")
else:
    sys.exit("Error: no 'design' found in config file. Check your configfile.\n")

#reference
if "reference" in config and config["reference"] in ["hg19", "hg38", "mm10", "mm9"]:
    sys.stderr.write("Reference set to: " + config["reference"] + ".\n")
else:
    sys.exit("Error in 'reference' set! Check your configfile. Reference can be: hg19, hg38, mm10 or mm9.\n")
    
#read design file
design=pd.read_table(config["design"], sep=",")

#check if input_format and variant_calling_mode are agree with design format
format_design=["sample_id","R1_fastq","R2_fastq"]
if set(format_design).issubset(design.columns):
  sys.stderr.write("Design file well formated.\n")
  design=design[format_design]
else: sys.exit("Error in the format of the design file: missing column(s).")

#check if all sample_id are differents
if not len(set(design["sample_id"])) == len(design["sample_id"]): sys.exit("Error: All sample_id have to be different! Check your design file.")

#check if all files are differents on each line
for i in range(0, len(design["sample_id"]), 1):
  if not len(set(design.iloc[i])) == len(design.iloc[i]): sys.exit("Error: All files have to be different! Check your design file.")

#make the data structure
SAMPLE_NAME = []
FQ_NAME = []
ORIG_FILE = []
SYMLINK_FILES = []
for col in design.drop(columns = "sample_id").columns.tolist():
  for line in range(0,len(design[col].tolist()),1):
    SAMPLE_NAME.append(design["sample_id"].iloc[line])
    ORIG_FILE.append(design[col].iloc[line])
    if col == "R1_fastq": SUPPL_NAME = "_R1"
    if col == "R2_fastq": SUPPL_NAME = "_R2"
    FQ_NAME.append(design["sample_id"].iloc[line] + SUPPL_NAME)
    SYMLINK_FILES.append(os.path.normpath(OUTPUT_DIR + "/symlink_input/" + design["sample_id"].iloc[line] + SUPPL_NAME + ".fastq.gz"))

#samples_order_for_ggplot
if not set(config["samples_order_for_ggplot"].split(",")).issubset(design["sample_id"]):
    sys.exit("Error sampes names of 'samples_order_for_ggplot' from config file don't match 'Sample_id' of design file.")



sys.stderr.write("\n########################### Run ############################\n\n")

### rule all ###

rule all:
    input:
        #symlink_qc_filtering
        expand(os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{sample_name}_R1.fastq.gz"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{sample_name}_R2.fastq.gz"), sample_name=SAMPLE_NAME),
        os.path.normpath(OUTPUT_DIR + "/multiqc.html"),
        #SPRINT
        expand(os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/{sample_name}.bam"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/PARAMETER.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/SPRINT_identified_all.res"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/SPRINT_identified_hyper.res"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/SPRINT_identified_regular.res"), sample_name=SAMPLE_NAME),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Number_of_all_edition_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Number_of_each_edition_by_sample_log_y_scale.png"),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Number_of_each_edition_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Number_of_reads_supporting_the_edition_for_each_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Percentages_of_reads_supporting_the_edition_for_each_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Summary_table_counts_SPRINT_by_strand.tsv"),
        os.path.normpath(OUTPUT_DIR + "/SPRINT/Summary_table_counts_SPRINT.tsv"),
        #RNAEditingIndexer
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/summary/EditingIndex.csv"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Coverage_by_base_reference_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Scores_of_all_edition_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Scores_of_edition_and_background_noise_for_each_edition_type.png")
    message:
        "Pipeline finished!"


### global wildcard constraints ###

wildcard_constraints:
    fq_name = '|'.join([x for x in FQ_NAME]),
    sample_name = '|'.join([x for x in SAMPLE_NAME])


### other rules ###

include: "rules/symlink_qc_filtering.smk"
include: "rules/SPRINT.smk"
include: "rules/RNAEditingIndexer.smk"







