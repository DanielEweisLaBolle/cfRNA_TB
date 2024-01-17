##-------------------------------
# PRE-PROCESS
# - filter low quality reads
# - remove duplicates either by sequence or UMI
# - remove three nucleotides from TSO

rule filter_q30:
    input:
        r1 = SAMPLE_LOC+'{sample}_R1.fastq.gz',
        r2 = SAMPLE_LOC+'{sample}_R2.fastq.gz'
    output:
        r1_clean = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.clean.fastq.gz',
        r2_clean = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.clean.fastq.gz'),
        r1_unpaired = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.unpaired.fastq.gz'),
        r2_unpaired = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.unpaired.fastq.gz')
    threads: trim_threads
    shell:
        """
        java -jar /programs/trimmomatic/trimmomatic-0.39.jar \
            PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.r1_clean} {output.r1_unpaired} {output.r2_clean} {output.r2_unpaired}  \
            AVGQUAL:30
        """ 
        
        
rule dedup:
    input:
        r1_clean = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.clean.fastq.gz',
        r2_clean = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.clean.fastq.gz'
    output:
        r1_dedup = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.dedup.fastq.gz',
        r2_dedup = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.dedup.fastq.gz')
    threads: 5
    params:
        lib_prep = lambda w: config['PREP'][w.sample],
        log_loc = OUTPUT+'{project}/sample_output/{sample}/{sample}.umitool.extract.log'
    shell:
        """
        /programs/bbmap-38.90/clumpify.sh in1={input.r1_clean} in2={input.r2_clean} out1={output.r1_dedup} out2={output.r2_dedup} dedupe
        """


rule headcrop:
    input: 
        r2_dedup = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.dedup.fastq.gz'
    output:
        r2_trim = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.headcrop.fastq.gz')
    threads: trim_threads
    params: 
        lib_prep = lambda w: config['PREP'][w.sample]
    shell:
        """
        java -jar /programs/trimmomatic/trimmomatic-0.39.jar \
            SE -threads {threads} \
            {input.r2_dedup} \
            {output.r2_trim}  \
            HEADCROP:3 
        """ 
        
rule bbmap_repair:
    input:
        r1_dedup = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.dedup.fastq.gz',
        r2_trim = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.headcrop.fastq.gz'
    output:
        r1_repaired = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.repaired.fastq.gz'),
        r2_repaired = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.repaired.fastq.gz')
    shell:
        """
        /programs/bbmap-38.90/repair.sh in1={input.r1_dedup} in2={input.r2_trim} out1={output.r1_repaired} out2={output.r2_repaired} repair
        """
##-------------------------------
# ALIGN AND COUNT
# - align using STAR with standard ENCODE parameters
# - count using RSEM
        
rule align:
    input:
        r1_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.repaired.fastq.gz',
        r2_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.repaired.fastq.gz'
    output:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
        bam_trscrpt = OUTPUT + '{project}/sample_output/{sample}/{sample}_Aligned.toTranscriptome.out.bam',
        STAR_log = OUTPUT+'{project}/sample_output/{sample}/{sample}_Log.final.out',
        u_1= OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R1.fastq',
        u_2 = OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R2.fastq',
    params:
        prefix = OUTPUT+"{project}/sample_output/{sample}/{sample}_",
        STAR = config["STAR"],
        STAR_IDX = lambda w: config["STAR_"+config['REF_DICT'][w.sample]],
        lib_prep = lambda w: config['PREP'][w.sample]
    threads: STAR_threads
    log: 'logs/{project}/{sample}.alignment.log'
    shell:
        """
            {params.STAR} \
            --genomeDir {params.STAR_IDX} \
            --readFilesIn {input.r1_repaired} {input.r2_repaired} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix}\
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --outSAMmultNmax 20 \
            --outSAMprimaryFlag AllBestScore \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignMatesGapMax 1000000 \
            --alignIntronMax 1000000 \
            --quantMode TranscriptomeSAM \
            --outReadsUnmapped Fastx 
            
            fastqc {input.r1_repaired}
            fastqc {input.r2_repaired}
            
            mv {params.prefix}Unmapped.out.mate1 {output.u_1}
            fastqc {output.u_1}
            mv {params.prefix}Unmapped.out.mate2 {output.u_2}
            fastqc {output.u_2}

        """

# --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \

            
rule RSEM:
    input:
        transcriptome_bam = OUTPUT + '{project}/sample_output/{sample}/{sample}_Aligned.toTranscriptome.out.bam'
    output:
        RSEM_bam_final = OUTPUT+'{project}/sample_output/{sample}/{sample}.RSEM.transcript.bam',
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}.RSEM.genes.results'
    params:
        RSEM_ref = lambda w: config["RSEM_"+config['REF_DICT'][w.sample] ],
        output_prefix = OUTPUT + '{project}/sample_output/{sample}/{sample}.RSEM'
    threads: RSEMCountThreads
    log: "logs/{project}/{sample}.featureCounts.log"
    params:
        lib_prep = lambda w: config['PREP'][w.sample]
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """        
        /programs/RSEM-1.3.3/bin/rsem-calculate-expression -p {threads} --strandedness reverse \
            --seed 42 --estimate-rspd --append-names --alignments --paired-end \
            {input.transcriptome_bam} {params.RSEM_ref} {params.output_prefix}
        """
        
        
##-------------------------------
# ALTERNATIVE COUNTING
        
rule feature_counts:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    output:
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.txt',
        summary = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.txt.summary'
    params:
        ANNOTATION = lambda w: config["ANNOTATION_"+config['REF_DICT'][w.sample] ]
    threads: featureCountThreads
    log: "logs/{project}/{sample}.featureCounts.log"
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """
            featureCounts \
            -p \
            -T {threads} \
            -s 2 \
            -O --fraction \
            -a {params.ANNOTATION} \
            -o {output.counts} \
            {input.bam}

            sed -i 's/\_sort\_dedup//g' {output.summary}
        """
        # -p: isPairedEnd 
        # -T: threads to use
        # -s: reverse stranded (2)
        # -O: allow multi overlap
        # --fraction: split read count between overlapping features

        
rule feature_counts_gene_name:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    output:
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.names.txt',
        summary = OUTPUT + '{project}/sample_output/{sample}/{sample}_feature_counts.names.txt.summary'
    params:
        ANNOTATION = lambda w: config["ANNOTATION_"+config['REF_DICT'][w.sample] ]
    threads: featureCountThreads
    log: "logs/{project}/{sample}.featureCounts.log"
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """
            featureCounts \
            -p \
            -T {threads} \
            -s 2 \
            -O --fraction \
            -a {params.ANNOTATION} \
            -g gene_name \
            -o {output.counts} \
            {input.bam}

            sed -i 's/\_sort\_dedup//g' {output.summary}
        """ 
        
rule htseq_count:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    output:
        counts = OUTPUT + '{project}/sample_output/{sample}/{sample}_htseq_counts.txt'
    log: 'logs/{project}/{sample}.htseq.log'
    params:
        ANNOTATION = lambda w: config["ANNOTATION_"+config['REF_DICT'][w.sample] ]
    conda:
        "../envs/cfRNA.yaml"
    threads: 10
    shell:
        """
        export PYTHONPATH=/programs/HTSeq-0.11.2/lib64/python3.6/site-packages/
        export PATH=/programs/HTSeq-0.11.2/bin:$PATH

        htseq-count \
                -f bam \
                -r pos \
                -s reverse \
                -m union \
                -i gene_name \
                --nonunique all \
                --secondary-alignments ignore \
                {input} \
                {params.ANNOTATION} > {output}      
        """
        
##-------------------------------
# ARCHIVED CODE
##-------------------------------


#         if [ {params.lib_prep} = "SMARTer_pico_v3" ]; then
            
#             samtools index {output.bam}
#             umi_tools dedup -I {output.bam} --output-stats={params.prefix}deduplicated --paired -S {output.bam}.tmp
#             mv {output.bam}.tmp {output.bam}
#             rm {output.bam}.bai
                
#             samtools sort {output.bam_trscrpt} > {output.bam_trscrpt}.sort
#             samtools index {output.bam_trscrpt}.sort
#             umi_tools dedup -I {output.bam_trscrpt}.sort --output-stats={params.prefix}.txn.deduplicated --paired --chimeric-pairs discard --unpaired-reads discard -S {output.bam_trscrpt}.dedup
#             /programs/RSEM-1.3.3/bin/convert-sam-for-rsem -p {threads} {output.bam_trscrpt}.dedup {output.bam_trscrpt}
                
#             rm {output.bam_trscrpt}.dedup
#             rm {output.bam_trscrpt}.sort
#             rm {output.bam_trscrpt}.sort.bai
            
            
#         fi
        
        
# --outSAMmultNmax 20 \                      <- int: max number of multiple alignments for a read that will be output to the SAM/BAM files
# --outSAMprimaryFlag AllBestScore \         <- which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG
# --outFilterMismatchNmax 999 \              <- maximum number of mismatches per pair, large number switches off this filter (999 turns it off)
# --outFilterMismatchNoverReadLmax 0.04 \    <- max number of mismatches per pair relative to read length: for 2x100b, max number of mis- matches is 0.04*200=8 for the paired read
# --alignIntronMin 20 \                      <- minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion
# --alignMatesGapMax 999                     <- maximum gap between two mates, if 0, max intron gap will be determined by
# --alignIntronMax 1000000 \                 <- maximum intron size, if 0, max intron size will be determined by
# --quantMode TranscriptomeSAM \
# --outReadsUnmapped Fastx \

# rule dedup_TMP:
#     input:
#         bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
#         bam_trscrpt = OUTPUT + '{project}/sample_output/{sample}/{sample}_Aligned.toTranscriptome.out.bam'
#     output:
#         dedup_txt = OUTPUT+'{project}/sample_output/{sample}/{sample}_dedup.txt'
#     params:
#         prefix = OUTPUT+"{project}/sample_output/{sample}/{sample}_",
#         lib_prep = lambda w: config['PREP'][w.sample]
#     threads: STAR_threads
#     log: 'logs/{project}/{sample}.alignment.log'
#     shell:
#         """
#         if [ {params.lib_prep} = "SMARTer_pico_v3" ]; then
            
                
#             samtools sort {input.bam_trscrpt} > {input.bam_trscrpt}.sort
#             samtools index {input.bam_trscrpt}.sort
#             umi_tools dedup -I {input.bam_trscrpt}.sort --output-stats={params.prefix}.txn.deduplicated --paired --chimeric-pairs discard --unpaired-reads discard -S {input.bam_trscrpt}.dedup
#             samtools sort -n {input.bam_trscrpt}.dedup > {input.bam_trscrpt}.dedup.nsort
#             /programs/RSEM-1.3.3/bin/convert-sam-for-rsem -p {threads} {input.bam_trscrpt}.dedup.nsort {input.bam_trscrpt}
                
#             rm {input.bam_trscrpt}.sort
#             rm {input.bam_trscrpt}.sort.bai
#             rm {input.bam_trscrpt}.dedup
#             rm {input.bam_trscrpt}.dedup.nsort
            
            
#         fi
        
#         touch {output.dedup_txt}
#         """