rule align_SARSCOV2:
    input:
        u_1= OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R1.fastq',
        u_2 = OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R2.fastq'
    output:
        covid_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_SARSCOV2_Aligned.sortedByCoord.out.bam'
    params:
        prefix = OUTPUT+"{project}/sample_output/{sample}/{sample}_SARSCOV2_",
        STAR = config["STAR"],
        STAR_IDX = '/workdir/cjl332/cfrna/references/covid19_fa/star'
    threads: STAR_threads
    log: 'logs/{project}/{sample}.SARSCOV2.alignment.log'
    shell:
        """
            {params.STAR} \
            --genomeDir {params.STAR_IDX} \
            --readFilesIn {input.u_1} {input.u_2} \
            --outFileNamePrefix {params.prefix}\
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --outSAMmultNmax 20 \
            --outSAMprimaryFlag AllBestScore \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 
        """