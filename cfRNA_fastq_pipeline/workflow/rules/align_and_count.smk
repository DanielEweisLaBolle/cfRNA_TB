##-------------------------------
# QUALITY TRIM
##-------------------------------
rule filter_q30:
    input:
        r1 = lambda w: get_fastq_path(w.sample) + '{sample}_R1.fastq.gz',
        r2 = lambda w: get_fastq_path(w.sample) + '{sample}_R2.fastq.gz'
    output:
        r1_clean = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.clean.fastq.gz',
        r2_clean = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.clean.fastq.gz')
    threads: TRIM_THREADS
    params:
        ref_apt = MPATH + config["TRUSEQ_ADAPT"],
        BBDUK = MPATH + config["BBDUK"],
        bbduk_output = OUTPUT+'{project}/sample_output/{sample}/{sample}_bbduk_check.out'
    shell:
        """
        {params.BBDUK} -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1_clean} out2={output.r2_clean} ref={params.ref_apt} ktrim=r k=21 mink=11 hdist=2 qtrim=rl trimq=10 tpe tbo >& {params.bbduk_output}
        """ 
        
##-------------------------------
# PREPARE READS TO BE ALIGNED
##-------------------------------
rule rename_headcrop_repair:
    input: 
        r1_clean = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.clean.fastq.gz',
        r2_clean = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.clean.fastq.gz'
    output:
        r1_named = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.named.fastq.gz'),
        r2_named = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.named.fastq.gz'),
        r2_trim = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.headcrop.fastq.gz'),
        r1_repaired = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.repaired.fastq.gz'),
        r2_repaired = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.repaired.fastq.gz')
    threads: TRIM_THREADS
    params: 
        lib_prep = lambda w: sample_information.loc[w.sample, "prep_type"],
        TRIMMOMATIC = config["TRIMMOMATIC"],
        BB_REPAIR = MPATH + config["BB_REPAIR"],
        log_loc = OUTPUT+'{project}/sample_output/{sample}/{sample}.umitool.extract.log'
    shell:
        """
        if [ {params.lib_prep} = "SMARTer_pico_v2" ]
        then
        
            cp {input.r1_clean} {output.r1_named}
            cp {input.r2_clean} {output.r2_named}
            
            java -jar {params.TRIMMOMATIC}\
                SE -threads {threads} \
                {output.r2_named} \
                {output.r2_trim}  \
                HEADCROP:3 
                
        fi

        if [ {params.lib_prep} = "SMARTer_pico_v3" ]
        then
        
            umi_tools extract --log {params.log_loc} --bc-pattern=NNNNNNNN \
            -I {input.r2_clean} -S {output.r2_named} \
            --read2-in={input.r1_clean} --read2-out={output.r1_named}

            java -jar {params.TRIMMOMATIC} \
                SE -threads {threads} \
                {output.r2_named} \
                {output.r2_trim}  \
                HEADCROP:6 
         
        fi
        
        {params.BB_REPAIR} in1={output.r1_named} in2={output.r2_trim} out1={output.r1_repaired} out2={output.r2_repaired} repair
        """ 
        
##-------------------------------
# ALIGN
##-------------------------------
rule align:
    input:
        r1_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.repaired.fastq.gz',
        r2_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.repaired.fastq.gz'
    output:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
        bam_trscrpt = OUTPUT + '{project}/sample_output/{sample}/{sample}_Aligned.toTranscriptome.out.bam',
        STAR_log = OUTPUT+'{project}/sample_output/{sample}/{sample}_Log.final.out',
        r1_unmapped = OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R1.fastq',
        r2_unmapped = OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R2.fastq'
    params:
        prefix = OUTPUT+"{project}/sample_output/{sample}/{sample}_",
        STAR = MPATH + config["STAR"],
#         STAR = "/programs/STAR-2.7.0f/bin/Linux_x86_64_static/STAR" , 
        STAR_IDX = MPATH + config["STAR_"+GENOME_REF],     # lambda w: config["STAR_"+config['REF_DICT'][w.sample]],
        lib_prep = lambda w: sample_information.loc[w.sample, "prep_type"]
    threads: ALIGN_THREADS
    log: 'logs/{project}/{sample}.alignment.log'
    shell:
        """
        {params.STAR} --runThreadN {threads} --genomeDir {params.STAR_IDX} --outSAMtype BAM SortedByCoordinate --readFilesIn {input.r1_repaired} {input.r2_repaired} --readFilesCommand zcat --outFileNamePrefix {params.prefix} --quantMode TranscriptomeSAM --outReadsUnmapped Fastx
        
        mv {params.prefix}Unmapped.out.mate1 {output.r1_unmapped}
        mv {params.prefix}Unmapped.out.mate2 {output.r2_unmapped}

        
        """
        
##-------------------------------
# DEDUPLICATE
##-------------------------------
rule dedup:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    output:
        dedup_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam',
        metrics = OUTPUT + '{project}/sample_output/{sample}/{sample}_dedup_metrics.txt'
    params:
        lib_prep = lambda w: sample_information.loc[w.sample, "prep_type"],
        prefix = OUTPUT+"{project}/sample_output/{sample}/{sample}_",
        PICARD = MPATH + config["PICARD"]
    threads: DEDUP_THREADS
    log: 'logs/{project}/{sample}.rmdup.log'
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """
        if [ {params.lib_prep} = "SMARTer_pico_v2" ]
        then
        
            java -XX:ParallelGCThreads={threads} -XX:ConcGCThreads={threads} -Xmx8000M -jar {params.PICARD} MarkDuplicates \
            I={input.bam} \
            O={output.dedup_bam} \
            M={output.metrics} \
            ASSUME_SORT_ORDER=coordinate \
            REMOVE_DUPLICATES=true
            
        fi

        if [ {params.lib_prep} = "SMARTer_pico_v3" ]
        then
        
            samtools index {input.bam}
            umi_tools dedup -I {input.bam} --output-stats={output.metrics} --paired -S {output.dedup_bam}
            rm {input.bam}.bai
            
            touch {output.metrics}
        
        fi
        """       

##-------------------------------
# COUNT FEATURES
##-------------------------------
rule feature_counts:
    input:
        bam =  OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam'
    output:
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.txt',
        summary = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.txt.summary',
        FTCOUNTS = OUTPUT + '{project}/sample_output/{sample}/{sample}_FTCOUNTS.txt'
    params:
        ANNOTATION = MPATH + config["ANNOTATION_"+GENOME_REF] # lambda w: config["ANNOTATION_"+config['REF_DICT'][w.sample] ]
    threads: FCOUNT_THREADS
    log: "logs/{project}/{sample}.featureCounts.DEDUPE.log"
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """
            featureCounts \
                --extraAttributes gene_name \
                -s 2 \
               -p -t exon -g gene_id \
               -T {threads} \
               -a {params.ANNOTATION} \
               -o {output.counts} \
                {input.bam}
                
            awk '{{ print $1,$8}}' {output.counts} > {output.FTCOUNTS}.tmp
            (echo "geneID {wildcards.sample} "; tail -n +3 {output.FTCOUNTS}.tmp) > {output.FTCOUNTS}
            rm {output.FTCOUNTS}.tmp
        """

##-------------------------------
# COUNT TARS
##-------------------------------

rule feature_counts_UTAR:
    input:
        bam =  OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam'
    output:
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.UTAR.txt',
        summary = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.UTAR.txt.summary'
    params:
        ANNOTATION = config["TAR_HG38_GTF"]
    threads: FCOUNT_THREADS
    log: "logs/{project}/{sample}.featureCounts.DEDUPE.log"
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """
            featureCounts \
                --extraAttributes gene_name \
                -s 2 \
               -p -t exon -g gene_id \
               -T {threads} \
               -a {params.ANNOTATION} \
               -o {output.counts} \
                {input.bam}
        """