rule phix_alignment:
    input:
        R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
        R2 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
    output:
        bam = temp('sample_output/phix/{sample}.bam'),
        unmapped_R1 = temp('sample_output/phix/unmapped/{sample}_unmapped_R1.fastq'),
        unmapped_R2 = temp('sample_output/phix/unmapped/{sample}_unmapped_R2.fastq')
    params:
        genome_prefix = 'references/phix/phix',
        unmapped_pe = 'sample_output/phix/unmapped/{sample}_unmapped_R%.fastq',
        prep_and_seq_type = get_sample_info
    shell:
        """
        bowtie2 -x {params.genome_prefix} \
                -1 {input.R1} -2 {input.R2} \
                --local --very-sensitive-local \
                --un-conc {params.unmapped_pe} | samtools view -bh - > {output.bam}
        """
