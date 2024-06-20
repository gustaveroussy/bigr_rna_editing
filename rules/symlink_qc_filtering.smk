"""
##########################################################################
These rules make the control-quality of fastq files
##########################################################################
"""

"""
This rule makes the symbolic links of fastq files with the good sample name.
"""

def symlink_rename_input_fq(wildcards):
    index = FQ_NAME.index(wildcards.fq_name)
    return ORIG_FILE[index]

rule symlink_rename_fq:
    input:
        fq = symlink_rename_input_fq
    output:
        fq_link = temp(os.path.normpath(OUTPUT_DIR + "/symlink_input/" + "/{fq_name}.fastq.gz"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 256),
        time_min = (lambda wildcards, attempt: attempt * 10)
    run:
        sys.stderr.write("\t Create symbolic link: \n")
        sys.stderr.write("\t From :" + "\t" + str(input.fq) + "\n")
        sys.stderr.write("\t To :" + "\t" + str(output.fq_link) + "\n")
        os.symlink(str(input.fq), str(output.fq_link))



"""
This rule makes the fastqc
"""

def get_fq(wildcards):
    index = FQ_NAME.index(wildcards.fq_name)
    return SYMLINK_FILES[index]

rule fastqc:
    input:
        fq_file = get_fq
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/fastqc/{fq_name}_fastqc.html")),
        temp(os.path.normpath(OUTPUT_DIR + "/fastqc/{fq_name}_fastqc.zip"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda :
        PIPELINE_DIR + "/envs/conda/fastQC.yaml"
    shell:
        """
        fastqc -o {OUTPUT_DIR}/fastqc/ {input.fq_file}
        """


"""
This rule makes the fastp
"""

def get_fq_R1(wildcards):
    index = FQ_NAME.index(wildcards.sample_name + "_R1")
    return SYMLINK_FILES[index]

def get_fq_R2(wildcards):
    index = FQ_NAME.index(wildcards.sample_name + "_R2")
    return SYMLINK_FILES[index]

rule fastp:
    input:
        R1_fq = get_fq_R1,
        R2_fq = get_fq_R2
    output:
        R1_fq = os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{sample_name}_R1.fastq"),
        R2_fq = os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{sample_name}_R2.fastq"),
        html = temp(os.path.normpath(OUTPUT_DIR + "/fastp/html/{sample_name}.fastp.html")),
        json = temp(os.path.normpath(OUTPUT_DIR + "/fastp/json/{sample_name}.fastp.json"))
    threads:
        5
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda :
        PIPELINE_DIR + "/envs/conda/fastp.yaml"
    shell:
        """
        #fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
        fastp --thread {threads} \
            --cut_front --cut_tail --cut_window_size 6 --cut_mean_quality 10 --unqualified_percent_limit 50 --n_base_limit 7 --average_qual 0 --length_required 15 --overrepresentation_analysis \
            --in1 {input.R1_fq} --in2 {input.R2_fq} \
            --out1 {output.R1_fq} --out2 {output.R2_fq} \
            --json {output.json} \
            --html {output.html}

        """


"""
This rule makes the multiqc to agregate all qc results into one html file.
"""

rule multiqc:
    input:
        expand(os.path.normpath(OUTPUT_DIR + "/fastqc/{fq_name}_fastqc.html"),fq_name=FQ_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/fastqc/{fq_name}_fastqc.zip"),fq_name=FQ_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/fastp/html/{sample_name}.fastp.html"),sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/fastp/json/{sample_name}.fastp.json"),sample_name=SAMPLE_NAME),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/fastqc/"))),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/fastp/html/"))),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/fastp/json/")))
    output:
        os.path.normpath(OUTPUT_DIR + "/multiqc_report.html"),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/multiqc_data")))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 5120),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda :
        PIPELINE_DIR + "/envs/conda/multiqc.yaml"
    shell:
        """
        cd {OUTPUT_DIR}
        multiqc ./fastp ./fastqc

        """
        