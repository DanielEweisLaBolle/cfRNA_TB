rule salmon:
    input:
        r1_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.repaired.fastq.gz',
        r2_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.repaired.fastq.gz'
    output:
        output = OUTPUT+'{project}/sample_output/{sample}/{sample}_salmon_quant.sf'
    params:
        prefix = OUTPUT+"{project}/sample_output/{sample}/{sample}_salmon_",
        SALMON_IDX = "/workdir/cjl332/cfrna/references/human/hg38/GRCh38.Salmon"
    threads: STAR_threads
    log: 'logs/{project}/{sample}.alignment.log'
    shell:
        """
        
        /programs/seqtk/seqtk shuffle -s 100 {input.r1_repaired} > {params.prefix}_RR1.fq.gz
        /programs/seqtk/seqtk shuffle -s 100 {input.r2_repaired} > {params.prefix}_RR2.fq.gz
        
        /programs/salmon-1.4.0/bin/salmon quant \
            --seqBias --gcBias --allowDovetail \
            --validateMappings \
            -p {threads} -l SR \
            -i {params.SALMON_IDX} \
            -1 {params.prefix}_RR1.fq.gz -2 {params.prefix}_RR2.fq.gz \
            --output {params.prefix}
            
        rm {params.prefix}_RR1.fq.gz
        rm {params.prefix}_RR2.fq.gz
        """