"""
##########################################################################
These rules make the editing analysis by RNAEditingIndexer
##########################################################################
"""

"""
This rule makes the RNAEditingIndexer analysis
"""
rule RNAEditingIndexer:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/{sample_name}/{sample_name}.sortedByCoord.bam")
    output:
        res = temp(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/tmp/{sample_name}/summary/EditingIndex.csv")),
        log = temp(directory(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/tmp/{sample_name}/log/")))
    threads:
        8
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 10240),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        bam_dir = (lambda wildcards: os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/input/"+ wildcards.sample_name)),
        cmpileup = (lambda wildcards: os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/tmp/"+ wildcards.sample_name + "/cmpileup/")),
        summary =(lambda wildcards: os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/tmp/"+ wildcards.sample_name + "/summary/")),
        ref = config["reference"],
        extra = config["RNAEditingIndexer_extra"] if "RNAEditingIndexer_extra" in config else ""
    shell:
        """
        singularity exec --no-home -B{OUTPUT_DIR} {PIPELINE_DIR}/envs/singularity/RNAEditingIndexer.simg \
        /bin/AEI/RNAEditingIndexer/RNAEditingIndex \
        {params.extra} \
        --ts 1 \
        --tsd {threads} \
        -d {params.bam_dir} \
        -f .sortedByCoord.bam \
        -l {output.log} \
        -o {params.cmpileup} \
        -os {params.summary} \
        --genome {params.ref} \
        --verbose && \
        rm -r {output.log}/flags && \
        mkdir -p {OUTPUT_DIR}"/RNAEditingIndexer/log/" && \
        for log_file in $(ls {output.log}); do mv {output.log}/$log_file {OUTPUT_DIR}/RNAEditingIndexer/log/; done

        """

"""
This rule merge the RNAEditingIndexer results
"""
rule RNAEditingIndexer_merge:
    input:
        expand(os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/tmp/{sample_name}/summary/EditingIndex.csv"),sample_name=list(dict.fromkeys(SAMPLE_NAME)))
    output:
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/summary/EditingIndex.csv")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 1024),
        time_min = (lambda wildcards, attempt: attempt * 60)
    shell:
        """
        head -n1 {input[0]} > {output}
        for result in {input}
        do
            tail -n +2 ${{result}} >> {output}
        done

        """


"""
This rule makes the RNAEditingIndexer graphs
"""
rule RNAEditingIndexer_summary:
    input:
        table = os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/summary/EditingIndex.csv")
    output:
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Coverage_by_base_reference_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Scores_of_all_edition_by_sample.png"),
        os.path.normpath(OUTPUT_DIR + "/RNAEditingIndexer/Scores_of_edition_and_background_noise_for_each_edition_type.png")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 1024),
        time_min = (lambda wildcards, attempt: attempt * 60)
    params:
        samples_order_for_ggplot = config["samples_order_for_ggplot"]
    shell:
        """
        singularity exec --no-home -B{PIPELINE_DIR},{OUTPUT_DIR} {PIPELINE_DIR}/envs/singularity/R_graphs.simg \
        Rscript {PIPELINE_DIR}/scripts/RNAEditingIndexer_summary_results.R --EditingIndex {input.table} --samples_order_for_ggplot {params.samples_order_for_ggplot}  

        """
