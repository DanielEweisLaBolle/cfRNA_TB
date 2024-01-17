

rule bbduk_filter:
    input:
        read1 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R1.fastq',
        read2 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R2.fastq',
        ref = '../references/hg38/hg38.fa'
    output:
        read1 = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1.fastq',
        read2 = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2.fastq'
    resources:
        mem_mb = 80000 #can change if necessary
    params:
        bbduk = config["BBDUK"]
    shell:
        """
          {params.bbduk} in1={input.read1} in2={input.read2} \
          out1={output.read1} out2={output.read2} \
          -Xmx{resources.mem_mb}m \
          prealloc=t \
          ref={input.ref} k=20
        """

rule trim_META:
    input: 
        r1 = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1.fastq',
        r2 = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1.fastq'
    output:
        r1_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1_trimmed.fq',
        r2_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2_trimmed.fq',
        r1_unpaired = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1_unpaired_trimmed.fq',
        r2_unpaired = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2_unpaired_trimmed.fq',
        trim_log = OUTPUT+'{project}/{sample}/{sample}_nonhuman_trimmomcatic.log'
    threads: trim_threads
    params:
        output_dir = OUTPUT+"{project}/{sample}/"
    log: 'logs/{project}/{sample}.trim.log'
    shell:
        """
        java -jar /programs/trimmomatic/trimmomatic-0.39.jar \
            PE -phred33 -threads {threads} -trimlog {output.trim_log}\
            {input.r1} {input.r2} \
            {output.r1_trim} {output.r1_unpaired} \
            {output.r2_trim} {output.r2_unpaired} \
            ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:36
        """        
		
rule filter_quality:
    input:
        r1_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1_trimmed.fq',
        r2_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2_trimmed.fq'
    output:
        r1_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1_filtered.fq',
        r2_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2_filtered.fq'
    params:
        SEQKIT = config["SEQKIT"]
    shell:
        """
        {params.SEQKIT} seq -Q 32 {input.r1_trim} > {output.r1_trim}
        {params.SEQKIT} seq -Q 32 {input.r2_trim} > {output.r2_trim}
        """

rule flash_merge:
    input:
        r1_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1_filtered.fq',
        r2_trim = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2_filtered.fq'
    output:
        o1 = OUTPUT + '{project}/{sample}/{sample}_tmp.merged.extendedFrags.fastq',
        o5 = temp(OUTPUT + '{project}/{sample}/{sample}_tmp.merged.notCombined_1.fastq'),
        o6 = temp(OUTPUT + '{project}/{sample}/{sample}_tmp.merged.notCombined_2.fastq'),
        o7 = temp(OUTPUT + '{project}/{sample}/{sample}.merged.fastq'),
        o8 = temp(OUTPUT + '{project}/{sample}/{sample}.R1.unmerged.fastq'),
        o9 = temp(OUTPUT + '{project}/{sample}/{sample}.R2.unmerged.fastq')
    threads:4
    params:
        dir = OUTPUT+'{project}/{sample}/',
        FLASH = config["FLASH"],
        prefix = '{sample}_'
    shell:
        """
        {params.FLASH} -m 10 {input.r1_trim} {input.r2_trim} -d {params.dir} -o {params.prefix}tmp.merged
        cat {output.o1} > {output.o7}
        cat {output.o5} > {output.o8}
        cat {output.o6} > {output.o9}
        """

rule filter_length:
    input:
        unpaired = OUTPUT + '{project}/{sample}/{sample}.merged.fastq',
        r1 = OUTPUT + '{project}/{sample}/{sample}.R1.unmerged.fastq',
        r2 = OUTPUT + '{project}/{sample}/{sample}.R2.unmerged.fastq'
    output:
        unpaired = OUTPUT + '{project}/{sample}/{sample}.merged.filtered.fastq',
        r1 = OUTPUT + '{project}/{sample}/{sample}.R1.filtered.fastq',
        r2 = OUTPUT + '{project}/{sample}/{sample}.R2.filtered.fastq'
    params:
        SEQKIT = config["SEQKIT"]
    shell:
        """
        {params.SEQKIT} seq -m 50 -M 140 {input.unpaired} > {output.unpaired}
        {params.SEQKIT} seq -m 50 -M 140 {input.r1} > {output.r1}
        {params.SEQKIT} seq -m 50 -M 140 {input.r2} > {output.r2}

        """
 
rule blast:
    input:
       unpaired = OUTPUT + '{project}/{sample}/{sample}.merged.filtered.fastq',
       r1 = OUTPUT + '{project}/{sample}/{sample}.R1.filtered.fastq',
       r2 = OUTPUT + '{project}/{sample}/{sample}.R2.filtered.fastq',
       db = '../references/stds/NCBIGenomes06.fna'
    output:
       bam1 = OUTPUT + '{project}/{sample}/{sample}.merged.blast',
       bam2 = OUTPUT + '{project}/{sample}/{sample}.R1.blast'
    threads: blast_threads
    params:
        HSBLASTN = config["HSBLASTN"]
    shell:
        """
        {params.HSBLASTN} align -query {input.unpaired} \
        -db {input.db} \
        -evalue 0.0001 \
        -perc_identity 95 \
        -num_threads {threads} \
        -outfmt 0 > {output.bam1}

        {params.HSBLASTN} align -query {input.r1} \
        -db {input.db} \
        -evalue 0.0001 \
        -perc_identity 95 \
        -num_threads {threads} \
        -outfmt 6 > {output.bam2}


        """
