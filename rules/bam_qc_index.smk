
"""
These rules sort and index bam files with samtools
"""
rule samtools_sort:
    input:
        os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/tmp/genome/all.bam")
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam"))
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    conda :
        PIPELINE_DIR + "/envs/conda/samtools.yaml"
    shell:
        """
        res_mem=$(({resources.mem_mb}-2000))
        echo $res_mem
        res_mem=$(echo $res_mem"M")
        samtools sort {input} -o {output} -@ {threads} -m $res_mem
        """

rule samtools_index:
    input:
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam")
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam.bai"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 250),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        PIPELINE_DIR + "/envs/conda/samtools.yaml"
    shell:
        "samtools index {input}"


"""
These rules compute qc on bam files with samtools
"""

rule samtools_stats:
    input:
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam")
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/Alignment_Quality_Control/samtools_stats/{sample_name}_stats.txt"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 250),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        PIPELINE_DIR + "/envs/conda/samtools.yaml"
    shell:
        "samtools stats {input} > {output}"

rule samtools_idxstats:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam"),
        bai = os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam.bai")
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/Alignment_Quality_Control/samtools_idxstats/{sample_name}_idxstats.txt"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 250),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        PIPELINE_DIR + "/envs/conda/samtools.yaml"
    shell:
        "samtools idxstats {input.bam} > {output}"

"""
This rule combines qc computed on bam files with multiqc
"""
rule bam_multiqc:
    input:
        expand(os.path.normpath(OUTPUT_DIR + "/Alignment_Quality_Control/samtools_stats/{sample_name}_stats.txt"), sample_name = SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Alignment_Quality_Control/samtools_idxstats/{sample_name}_idxstats.txt"), sample_name = SAMPLE_NAME)
    output:
        os.path.normpath(OUTPUT_DIR + "/multiqc_alignment_report.html"),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/multiqc_alignment_report_data")))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 250),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        PIPELINE_DIR + "/envs/conda/multiqc.yaml"
    shell:
        """
        cd {OUTPUT_DIR}/ &&
        multiqc ./Alignment_Quality_Control/ --filename multiqc_alignment_report.html -c {PIPELINE_DIR}/config/bam_multiqc_config.yaml
        """
