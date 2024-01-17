rule BP_decon:
    input:
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.txt',
    output:
        celltype_frac = OUTPUT + '{project}/sample_output/{sample}/{sample}_BP_celltypes.tsv',
        celltype_frac_pc = OUTPUT + '{project}/sample_output/{sample}/{sample}_BP_celltypes.protein_coding.tsv'
    params:
        decon_script = config["BP_DECON"],
        biotype_key = MPATH + config["BIOTYPE_KEY_"+GENOME_REF],
        decon_ref = config['DECON_REF_COUNTS'],
        decon_meta = config['DECON_REF_META']
    threads: DECON_THREADS
    shell:
        """
        Rscript {params.decon_script} {input.counts} {params.decon_ref} {params.decon_meta} {params.biotype_key} {output.celltype_frac} {output.celltype_frac_pc} {threads} {wildcards.sample}
        """
    

    