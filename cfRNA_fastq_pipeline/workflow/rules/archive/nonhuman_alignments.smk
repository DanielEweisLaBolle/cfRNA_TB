rule align_covid:
    input:
        u_1= OUTPUT+'{project}/{sample}/{sample}_unmapped_R1.fastq',
        u_2 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R2.fastq'
    output:
        covid_bam = OUTPUT + '{project}/{sample}/{sample}_covid.bam',
        covid_btw_report = OUTPUT + '{project}/{sample}/{sample}_covid.bt2.txt'
    threads: alignment_thread
    params:
        covid_REF = config["COVID_REF"]
    shell:
        """
        bwa-mem -x {params.covid_REF} \
            -1 {input.u_1} -2 {input.u_2} \
            -p {threads} -q \
            --dovetail 2> {output.covid_btw_report} | samtools view -f3 -bh - > {output.covid_bam}
        """
		
# 		        bowtie2 -x {params.covid_REF} \
#             -1 {input.u_1} -2 {input.u_2} \
#             -p {threads} -q \
#             --dovetail 2> {output.covid_btw_report} | samtools view -f3 -bh - > {output.covid_bam}