rule filter_contam:
    input:
        tblatpe = 'sample_output/blast/consolidated_blast/{sample}.tblat.pe',
        reference_fasta = 'databases/GenomeDB/NCBIGenomes06.fna',
        LUT = 'LUTGrammy/taxids_names_lengths_tax.tab',
        unmapped_R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
        unmapped_R2 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz',
        model1= 'CFS/CFSBS_MODEL_R1.joblib',
        model2= 'CFS/CFSBS_MODEL_R2.joblib'
    output:
        names = temp('sample_output/blast/consolidated_blast/{sample}.reads'),
        R1_sort = temp('sample_output/blast/consolidated_blast/{sample}_temp_sort_R1.fastq'),
        R2_sort = temp('sample_output/blast/consolidated_blast/{sample}_temp_sort_R2.fastq'),
        tblatpe_sort = temp('sample_output/blast/consolidated_blast/{sample}.tblat.sort.pe'),
        removed = 'sample_output/filtered_grammy/{sample}/removed.{sample}.blast',
        kept = 'sample_output/filtered_grammy/{sample}/{sample}.blast',
        tblat = 'sample_output/filtered_grammy/{sample}/{sample}.tblat.1'
    params:
        dir = 'sample_output/filtered_grammy/{sample}'
    log:
        'logs/contamination_filtering/{sample}.log'
    threads: 1
    shell:
        """
        cut -f2 {input.tblatpe} | tail -n +2 | cut -f2 -d'_' > {output.names}

        zcat {input.unmapped_R1} | cut -f1 -d'_' | cut -f1 -d' ' | \
            ./software/seqtk/seqtk subseq - {output.names} | \
            paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > {output.R1_sort}

        zcat {input.unmapped_R2} | cut -f1 -d'_' | cut -f1 -d' ' | \
            ./software/seqtk/seqtk subseq - {output.names} | \
            paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > {output.R2_sort}

        (head -n1 {input.tblatpe} && tail -n +2 {input.tblatpe} | sort -k2) > {output.tblatpe_sort}

        python ./scripts/contamination_filtering/filter3_004.py \
            -fn {output.R1_sort} {output.R2_sort} {output.tblatpe_sort} \
            -ft tblat \
            -mo {input.model1} {input.model2} \
            -cf f \
            -rt r \
            -d {params.dir} \
            -r {log} \
            -o {wildcards.sample}.blast \
            -mr {input.reference_fasta} \
            -LUT {input.LUT} -j {threads}
        tail -n +2 {output.kept} > {output.tblat}
        """
