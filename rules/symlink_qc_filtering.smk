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

rule raw_fastqc:
    input:
        fq_file = get_fq
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/raw_fastqc/{fq_name}_fastqc.html")),
        temp(os.path.normpath(OUTPUT_DIR + "/raw_fastqc/{fq_name}_fastqc.zip"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda :
        PIPELINE_DIR + "/envs/conda/fastqc.yaml"
    shell:
        """
        fastqc -o {OUTPUT_DIR}/raw_fastqc/ {input.fq_file}
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
        R1_fq = temp(os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{sample_name}_R1.fastq")),
        R2_fq = temp(os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{sample_name}_R2.fastq")),
        html = temp(os.path.normpath(OUTPUT_DIR + "/fastp/html/{sample_name}.fastp.html")),
        json = temp(os.path.normpath(OUTPUT_DIR + "/fastp/json/{sample_name}.fastp.json"))
    threads:
        6
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 1024),
        time_min = (lambda wildcards, attempt: attempt * 240)
    conda :
        PIPELINE_DIR + "/envs/conda/fastp.yaml"
    shell:
        """
        fastp --thread {threads} \
            --cut_front --cut_tail --cut_window_size 6 --cut_mean_quality 10 --unqualified_percent_limit 50 --n_base_limit 7 --average_qual 0 --length_required 15 --overrepresentation_analysis \
            --in1 {input.R1_fq} --in2 {input.R2_fq} \
            --out1 {output.R1_fq} --out2 {output.R2_fq} \
            --json {output.json} \
            --html {output.html}
        """

if SAMPLING_BOOL :
    """
    This rule makes the sampling of reads.
    """
    def get_input_seqtk(wildcards):
        if config["nb_sampled_reads"] == "auto":
            return expand(os.path.normpath(OUTPUT_DIR + "/fastp/json/{sample_name}.fastp.json"), sample_name = SAMPLE_NAME)
        else:
            return []
    def get_nb_sampled_reads(wildcards):
        if config["nb_sampled_reads"] == "auto":
            all_nb_reads_after_filtering=[]
            for file in get_input_seqtk(wildcards):
                with open(file) as f:
                    lines = f.readlines()
                    all_nb_reads_after_filtering.append(int(lines[16].replace('\t\t\t"total_reads":', '').replace(',\n', '')))
            config["nb_sampled_reads"] = max(min(all_nb_reads_after_filtering),config["min_nb_sampled_reads_for_auto"])
        return config["nb_sampled_reads"]/2
    rule seqtk:
        input:
            fq = os.path.normpath(OUTPUT_DIR + "/fastp/trimmed/{fq_name}.fastq"),
            all_json = get_input_seqtk
        output:
            temp(os.path.normpath(OUTPUT_DIR + "/fastp/trimmed_sampled/{fq_name}.fastq"))
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: attempt * 40960),
            time_min = (lambda wildcards, attempt: attempt * 720)
        params:
            nb_sampled_reads = get_nb_sampled_reads
        conda :
            PIPELINE_DIR + "/envs/conda/seqtk.yaml"
        shell:
            """
            seqtk sample -s100 {input.fq} {params.nb_sampled_reads} > {output}
            """
    """
    This rule makes the fastqc of sampled reads
    """
    rule sampled_fastqc:
        input:
            os.path.normpath(OUTPUT_DIR + "/fastp/trimmed_sampled/{fq_name}.fastq")
        output:
            temp(os.path.normpath(OUTPUT_DIR + "/sampled_fastqc/{fq_name}_fastqc.html")),
            temp(os.path.normpath(OUTPUT_DIR + "/sampled_fastqc/{fq_name}_fastqc.zip"))
        threads:
            1
        resources:
            mem_mb = (lambda wildcards, attempt: attempt * 4096),
            time_min = (lambda wildcards, attempt: attempt * 720)
        conda :
            PIPELINE_DIR + "/envs/conda/fastqc.yaml"
        shell:
            """
            fastqc -o {OUTPUT_DIR}/sampled_fastqc/ {input}
            """


"""
This rule makes the multiqc to agregate all qc results into one html file.
"""
def get_fq_multiqc_input(wildcards):
    if SAMPLING_BOOL :
       return expand(os.path.normpath(OUTPUT_DIR + "/sampled_fastqc/{fq_name}_fastqc.{ext}"),fq_name=FQ_NAME,ext=["html","zip"])
    else:
        return []
def get_multiqc_config(wildcards):
    if SAMPLING_BOOL:
        return "-c " + PIPELINE_DIR + "/config/fq_multiqc_sampled_config.yaml"
    else:
        return "-c " + PIPELINE_DIR + "/config/fq_multiqc_config.yaml"
def get_multiqc_path(wildcards):
    if SAMPLING_BOOL:
        return "./fastp ./raw_fastqc ./sampled_fastqc"
    else:
        return "./fastp ./raw_fastqc"
rule fq_multiqc:
    input:
        expand(os.path.normpath(OUTPUT_DIR + "/raw_fastqc/{fq_name}_fastqc.html"),fq_name=FQ_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/raw_fastqc/{fq_name}_fastqc.zip"),fq_name=FQ_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/fastp/html/{sample_name}.fastp.html"),sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/fastp/json/{sample_name}.fastp.json"),sample_name=SAMPLE_NAME),
        get_fq_multiqc_input
    output:
        os.path.normpath(OUTPUT_DIR + "/multiqc_raw_data_report.html"),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/multiqc_raw_data_report_data")))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 250),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        multiqc_config = get_multiqc_config,
        multiqc_path = get_multiqc_path
    conda :
        PIPELINE_DIR + "/envs/conda/multiqc.yaml"
    shell:
        """
        cd {OUTPUT_DIR}
        multiqc {params.multiqc_path} --filename multiqc_raw_data_report.html {params.multiqc_config}
        
        """
        
