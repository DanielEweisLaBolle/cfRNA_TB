# rule references_for_decontamination: #rule is necessary because bbduk only looks at either Forward strand or F and R (not just Reverse)
#     input:
#         hg19_CT = ancient('references/hg19_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa'),
#         hg19_GA = ancient('references/hg19_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'),
#         mic_ctl_CT = ancient('references/controls_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa'),
#         mic_ctl_GA = ancient('references/controls_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'),
#         phix_CT = ancient('references/phix_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa'),
#         phix_GA = ancient('references/phix_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'),
#         phix = ancient('references/phix/phix.fa'),
#         mic = ancient('references/controls_Bismark/mic.fa')
#     output:
#         CT = 'references/for_decontamination/CT.fa',
#         GA = 'references/for_decontamination/GA.fa',
#         both = 'references/for_decontamination/both.fa'
#     shell:
#         """
#         cat {input.phix} {input.mic} > {output.both}
#         cat {input.hg19_CT} {input.mic_ctl_CT} {input.phix_CT} > {output.CT}
#         cat {input.hg19_GA} {input.mic_ctl_GA} {input.phix_GA} | {SEQTK} seq -r - > {output.GA}
#         """

rule decontaminate_wgbs:
    input:
        unmapped_R1 = 'sample_output/phix/unmapped/{sample}_unmapped_R1.fastq',
        unmapped_R2 = 'sample_output/phix/unmapped/{sample}_unmapped_R2.fastq',
        CT = 'references/for_decontamination/CT.fa',
        GA = 'references/for_decontamination/GA.fa',
        both = 'references/for_decontamination/both.fa',
        mic_ctl_CT = 'references/controls_Bismark/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa',
        mic_ctl_GA = 'references/controls_Bismark/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa'
    output:
        r1_CT = temp('sample_output/decontaminate_wgbs/{sample}.CT_R1.fastq'),
        r2_CT = temp('sample_output/decontaminate_wgbs/{sample}.CT_R2.fastq'),
        r1_pass1 = temp('sample_output/decontaminate_wgbs/{sample}.pass1_R1.fastq'),
        r2_pass1 = temp('sample_output/decontaminate_wgbs/{sample}.pass1_R2.fastq'),
        r1_pass2 = temp('sample_output/decontaminate_wgbs/{sample}.pass2_R1.fastq'),
        r2_pass2 = temp('sample_output/decontaminate_wgbs/{sample}.pass2_R2.fastq'),
        decon_r1 = temp('sample_output/decontaminate_wgbs/{sample}_decontaminated_R1.fastq'),
        decon_r2 = temp('sample_output/decontaminate_wgbs/{sample}_decontaminated_R2.fastq'),
        nonhumanfa_r1 = temp('sample_output/nonhuman_fasta_wgbs/R1/{sample}_R1.fa'),
        nonhumanfa_r2 = temp('sample_output/nonhuman_fasta_wgbs/R2/{sample}_R2.fa')
    resources:
        mem_mb=80000
    shell:
        """
        python scripts/metagenome/insilico_conversion2.py -i {input.unmapped_R1} -o {output.r1_CT} -r R1
        python scripts/metagenome/insilico_conversion2.py -i {input.unmapped_R2} -o {output.r2_CT} -r R2
        {BBDUK} in1={output.r1_CT} in2={output.r2_CT} \
                out1={output.r1_pass1} out2={output.r2_pass1} \
                -Xmx{resources.mem_mb}m \
                prealloc=t rcomp=f \
                ref={input.GA} k=50
        {BBDUK} in1={output.r1_pass1} in2={output.r2_pass1} \
                out1={output.r1_pass2} out2={output.r2_pass2} \
                -Xmx{resources.mem_mb}m \
                prealloc=t rcomp=f \
                ref={input.CT} k=50
        {BBDUK} in1={output.r1_pass2} in2={output.r2_pass2} \
                out1={output.decon_r1} out2={output.decon_r2} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={input.both} k=20
        fastq_to_fasta -Q33 -i {output.decon_r1} -o {output.nonhumanfa_r1}
        fastq_to_fasta -Q33 -i {output.decon_r2} -o {output.nonhumanfa_r2}
        """

rule decontaminate_stds:
    input:
        unmapped_R1 = 'sample_output/phix/unmapped/{sample}_unmapped_R1.fastq',
        unmapped_R2 = 'sample_output/phix/unmapped/{sample}_unmapped_R2.fastq',
        genome = get_reference_genome_fasta,
        both = 'references/for_decontamination/both.fa'
    output:
        r1_pass1 = temp('sample_output/decontaminate_stds/{sample}.pass1_R1.fastq'),
        r2_pass1 = temp('sample_output/decontaminate_stds/{sample}.pass1_R2.fastq'),
        decon_r1 = temp('sample_output/decontaminate_stds/{sample}_decontaminated_R1.fastq'),
        decon_r2 = temp('sample_output/decontaminate_stds/{sample}_decontaminated_R2.fastq'),
        nonhumanfa_r1 = temp('sample_output/nonhuman_fasta_stds/R1/{sample}_R1.fa'),
        nonhumanfa_r2 = temp('sample_output/nonhuman_fasta_stds/R2/{sample}_R2.fa')
    resources:
        mem_mb=80000
    shell:
        """
        {BBDUK} in1={input.unmapped_R1} in2={input.unmapped_R2} \
                out1={output.r1_pass1} out2={output.r2_pass1} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={input.genome} k=50
        {BBDUK} in1={output.r1_pass1} in2={output.r2_pass1} \
                out1={output.decon_r1} out2={output.decon_r2} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={input.both} k=20
        fastq_to_fasta -Q33 -i {output.decon_r1} -o {output.nonhumanfa_r1}
        fastq_to_fasta -Q33 -i {output.decon_r2} | {SEQTK} seq -r - > {output.nonhumanfa_r2}
        """
