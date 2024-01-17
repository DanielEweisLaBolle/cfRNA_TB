rule CT_mismatch:
    input:
        tblatpe= 'sample_output/blast/consolidated_blast/{sample}.tblat.pe',
        reference_fasta = 'databases/GenomeDB/NCBIGenomes06.fna',
        unmapped_R1 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz',
        unmapped_R2 = 'sample_output/aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
    output:
        names = temp('sample_output/blast/consolidated_blast/{sample}.reads'),
        R1_sort = temp('sample_output/blast/consolidated_blast/{sample}_temp_sort_R1.fastq'),
        R2_sort = temp('sample_output/blast/consolidated_blast/{sample}_temp_sort_R2.fastq'),
        tblatpe_sort = temp('sample_output/blast/consolidated_blast/{sample}.tblat.sort.pe'),
        tblat1_CT = 'sample_output/filter_unconverted/grammy/{sample}.tblat.1'
    shell:
        """
        cut -f2 {input.tblatpe} | tail -n +2 | cut -f1 -d'_' > {output.names}

        zcat {input.unmapped_R1} | cut -f1 -d'_' | \
            ./software/seqtk/seqtk subseq - {output.names} | \
            paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > {output.R1_sort}

        zcat {input.unmapped_R2} | cut -f1 -d'_' | \
            ./software/seqtk/seqtk subseq - {output.names} | \
            paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > {output.R2_sort}


        tail -n +2 {input.tblatpe} | sort -k2 > {output.tblatpe_sort}
        python scripts/filter_conversion/filter2.py \
                                                --filenames {output.R1_sort} {output.R2_sort} {output.tblatpe_sort} \
                                                -ft n -c 1 \
                                                -rt t \
                                                --output_tblat {output.tblat1_CT} \
                                                --metagenomic_reference {input.reference_fasta}
        """

rule touchup:
    input:
        tblat1_CT = 'sample_output/filter_unconverted/grammy/{sample}.tblat.1',
        LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
    output:
        out = 'sample_output/filter_unconverted/grammy/{sample}.filtered.tmp'
    shell:
        """
        Rscript ./scripts/metagenome/contaminant_finder.R {input.tblat1_CT} {output.out} {input.LUT}
        """

#rule creat_ML_classi:
#    input:
#        'CFS/training_samples.tsv',
#        expand('sample_output/filter_unconverted/grammy/{training_sample}.filtered.tmp', training_sample = aget_training_samples)
#    output:
#        'CFS/model.pkl',
#        'CFS/model_roc.png'
#    threads:
#        5
#    shell:
#        """
#        python ./scripts/filter_conversion/create_model2.py --input {input[0]} --output {output[0]} -p {threads} --figure {output[1]}
#        """

rule contam_finder:
    input:
        filtered = 'sample_output/filter_unconverted/grammy/{sample}.filtered.tmp',
        model = 'CFS/model.pkl'
    output:
        filtered = 'sample_output/filter_unconverted/grammy/{sample}/{sample}.filtered',
        filt_tblat = 'sample_output/filter_unconverted/grammy/{sample}/{sample}.tblat.1'
    shell:
        """
        python ./scripts/filter_conversion/make_predictions.py --input_data {input.filtered} --model {input.model} --output_data {output.filtered} --output_grammy {output.filt_tblat}
        """
