rule runKraken2:
    input:
        r1_unmapped = OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R1.fastq',
        r2_unmapped = OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_R2.fastq'
    output: 
        r1_dedup = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_DEDUPE_R1.fastq'),
        r2_dedup = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_unmapped_DEDUPE_R2.fastq'),
        k_output = OUTPUT+'{project}/sample_output/{sample}/{sample}_kraken.out',
        k_report = OUTPUT+'{project}/sample_output/{sample}/{sample}_kraken.report'
    threads: 5
    params:
        CLUMPIFY = MPATH + config["CLUMPIFY"],
        kraken_db = MPATH + config["KRAKEN_DB"]
    shell:
        """
        {params.CLUMPIFY} in1={input.r1_unmapped} in2={input.r2_unmapped} out1={output.r1_dedup} out2={output.r2_dedup} dedupe
        
        kraken2 --db {params.kraken_db} \
        --threads {threads} --paired --use-names \
        --report {output.k_report} \
        {output.r1_dedup} {output.r2_dedup} > {output.k_output}
        """

rule bracken:
    input:
        kraken_report=OUTPUT+'{project}/sample_output/{sample}/{sample}_kraken.report'
    output:
        bracken_out=OUTPUT+"{project}/sample_output/{sample}/{sample}.bracken"
    params:
        kraken_db = MPATH + config["KRAKEN_DB"]
    shell:
        """
        bracken -d {params.kraken_db} -i {input.kraken_report} -o {output.bracken_out} -r 75
        """
  
rule extend_bracken:
    input:
        bracken_out=OUTPUT+"{project}/sample_output/{sample}/{sample}.bracken"
    output:
        bracken_extended=OUTPUT+"{project}/sample_output/{sample}/{sample}.bracken.extended"
    params:
        EX_BRACKEN = config["EX_BRACKEN"]
    shell:
        """
        {params.EX_BRACKEN} --bracken_file {input.bracken_out} --output {output.bracken_extended} 
        """