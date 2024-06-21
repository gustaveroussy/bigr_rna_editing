"""
##########################################################################
These rules make the editing analysis by RNAEditingIndexer
##########################################################################
"""


"""
This rule makes the sorting of bam files with samtools
"""
rule samtools_sort:
    input:
        bam=os.path.normpath(OUTPUT_DIR + "/SPRINT/{sample_name}/tmp/genome/all.bam")
    output:
        bam=temp(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}.sortedByCoord.bam"))
    threads:
        8
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda :
        PIPELINE_DIR + "/envs/conda/samtools.yaml"
    shell:
        """
        samtools sort {input.bam} -o {output.bam} -@ {threads} -m {resources.mem_mb}M
        """


"""
This rule makes the RNAEditingIndexer analysis
"""
rule RNAEditingIndexer:
    input:
        bam_list=expand(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}.sortedByCoord.bam"),sample_name=SAMPLE_NAME)
    output:
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/summary/EditingIndex.csv")
    threads:
        16
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        bam_dir=os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/"),
        log=os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/log/"),
        cmpileup=os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/cmpileup/"),
        summary=os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/summary/"),
        ref=config["reference"],
        extra=config["RNAEditingIndexer_extra"] if "RNAEditingIndexer_extra" in config else ""
    shell:
        """
        singularity exec --no-home -B{OUTPUT_DIR} {PIPELINE_DIR}/envs/singularity/RNAEditingIndexer.simg \
        /bin/AEI/RNAEditingIndexer/RNAEditingIndex \
        {params.extra} \
        --ts 1 \
        --tsd {threads} \
        -d {params.bam_dir} \
        -f .sortedByCoord.bam \
        -l {params.log} \
        -o {params.cmpileup} \
        -os {params.summary} \
        --genome {params.ref} \
        --verbose || rm -r {OUTPUT_DIR}/RNAEditingIndexer/output

        """


"""
This rule makes the RNAEditingIndexer graphs
"""
rule RNAEditingIndexer_summary:
    input:
        table=os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/summary/EditingIndex.csv")
    output:
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Coverage_by_base_reference_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Scores_of_all_edition_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Scores_of_edition_and_background_noise_for_each_edition_type.png")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 2048),
        time_min = (lambda wildcards, attempt: attempt * 60)
    params:
        samples_order_for_ggplot=config["samples_order_for_ggplot"]
    shell:
        """
        singularity exec --no-home -B{PIPELINE_DIR},{OUTPUT_DIR} {PIPELINE_DIR}/envs/singularity/R_graphs.simg \
        Rscript {PIPELINE_DIR}/scripts/RNAEditingIndexer_summary_results.R --EditingIndex {input.table} --samples_order_for_ggplot {params.samples_order_for_ggplot}  

        """
