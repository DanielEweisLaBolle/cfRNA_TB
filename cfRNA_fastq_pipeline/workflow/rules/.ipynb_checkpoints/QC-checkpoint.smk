##-------------------------------
# FASTQC
##-------------------------------
rule fastqc:
    input:
        r1 = lambda w: get_fastq_path(w.sample) + '{sample}_R1.fastq.gz',
        r2 = lambda w: get_fastq_path(w.sample) + '{sample}_R2.fastq.gz'
    output:
        OUTPUT + "{project}/sample_output/{sample}/{sample}_R1_fastqc.html",
        OUTPUT + "{project}/sample_output/{sample}/{sample}_R2_fastqc.html"
    threads: QC_THREADS
    params:
        outdir = OUTPUT + "{project}/sample_output/{sample}/"
    log: "logs/{project}fastqc/{sample}.fastqc.log"
    shell:
        """
            fastqc {input.r1} {input.r2} -t {threads} --outdir {params.outdir} -q
        """
            
##-------------------------------
# QUALIMAP
##-------------------------------
rule qualimap:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam',
    output:
        output_qualimap = OUTPUT + "{project}/sample_output/{sample}/rnaseq_qc_results.txt"
    params:
        outdir = OUTPUT+'{project}/sample_output/{sample}/',
        ANNOTATION = MPATH + config["ANNOTATION_" + GENOME_REF],
        QUALIMAP =  config["QUALIMAP"]
    threads: QC_THREADS
    shell:
        """
        {params.QUALIMAP} rnaseq \
            -bam {input} \
            -gtf {params.ANNOTATION} \
            --outdir {params.outdir} \
            --outfile {wildcards.sample}_qualimap \
            --sequencing-protocol strand-specific-reverse \
            --paired \
            --java-mem-size=10G
        """

        
##-------------------------------
# DNA CONTAMINATION
##-------------------------------
rule intron_exon_bams:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam'
    output:
        stats_file = OUTPUT+'{project}/sample_output/{sample}/{sample}_exon_intron.stats',
        exon_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_exon.bam',
        non_exon_bam = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_NONexon.bam'),
        intron_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_intron.bam',
        non_intron_exon_bam = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_NONexonintron.bam'),
        intergenic_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_intergenic.bam'
    params:
        exon_bed = MPATH + config["EXONS_" + GENOME_REF], # lambda w: config["EXONS_"+config['REF_DICT'][w.sample] ],
        intron_bed = MPATH + config["INTRONS_" + GENOME_REF], # lambda w: config["INTRONS_"+config['REF_DICT'][w.sample] ],
        intergenic_bed = MPATH + config["INTERGENIC_" + GENOME_REF] # lambda w: config["INTERGENIC_"+config['REF_DICT'][w.sample] ]
    threads: QC_THREADS
    shell:
        """
        samtools view -bf 0x2 -F 0x100 -U {output.non_exon_bam} -L {params.exon_bed} {input.bam} > {output.exon_bam}
        samtools view -c {output.exon_bam} >> {output.stats_file}
        
        samtools view -bf 0x2 -F 0x100 -U {output.non_intron_exon_bam} -L {params.intron_bed} {output.non_exon_bam} > {output.intron_bam}
        samtools view -c {output.intron_bam} >> {output.stats_file}
        
        samtools view -bf 0x2 -F 0x100 -L {params.intergenic_bed} {output.non_intron_exon_bam} > {output.intergenic_bam}
        samtools view -c {output.intergenic_bam} >> {output.stats_file}
        """      

##-------------------------------
# SPLIT BAMS
##-------------------------------
rule split_bams:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam'
    output:
        nuclear_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.NUCLEAR.bam',
        mito_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.MT.bam',
        ribo_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.RIBO.bam'
    params:
        rRNA_bed = MPATH + config["rRNA_HG38_BED"],
        ORGANISM = GENOME_REF # lambda w: config['REF_DICT'][w.sample]
    threads: QC_THREADS
    shell:
        """
        samtools index -@ {threads} {input.bam}
        if [[ {params.ORGANISM} = "GRCm38" ]]
        then
            samtools view -@ {threads} -F 0x100 -b -o {output.nuclear_bam} {input.bam} `seq 1 19 | sed 's/^/GRCm38_/'` 
            samtools view -@ {threads} -F 0x100 -b -o {output.mito_bam} {input.bam} GRCm38_MT 
        else
            samtools view -@ {threads} -F 0x100 -b -o {output.nuclear_bam} {input.bam} `seq 1 22 | sed 's/^/chr/'`
            samtools view -@ {threads} -F 0x100 -b -o {output.mito_bam} {input.bam} chrM 
        fi
        samtools view -@ {threads} -F 0x100 -b -o {output.ribo_bam} {input.bam} -L {params.rRNA_bed}
        
        """
        
##-------------------------------
# COUNT BIOTYPES
##-------------------------------
rule biotype_counts:
    input:
        bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam'
    output:
        counts = temp(OUTPUT + '{project}/sample_output/{sample}/{sample}_BIOTYPE_counts.txt'),
        summary = temp(OUTPUT + '{project}/sample_output/{sample}/{sample}_BIOTYPE_counts.txt.summary'),
        BIOTYPE = OUTPUT + '{project}/sample_output/{sample}/{sample}_BIOTYPES.txt'
    params:
        ANNOTATION = MPATH + config["ANNOTATION_AUTO_" + GENOME_REF] # lambda w: config["ANNOTATION_AUTO_"+config['REF_DICT'][w.sample] ]
    threads: FCOUNT_THREADS
    log: "logs/{project}/{sample}.featureCounts.log"
    conda:
        "../envs/cfRNA.yaml"
    shell:
        """
        featureCounts \
            -s 2 \
            -p -t exon -g gene_type \
            -T {threads} \
            -a {params.ANNOTATION} \
            -o {output.counts} \
             {input.bam}
                  
        awk '{{ print $1,$7}}' {output.counts} > {output.BIOTYPE}.tmp
        (echo "biotype {wildcards.sample}"; tail -n +3 {output.BIOTYPE}.tmp) > {output.BIOTYPE}
        rm {output.BIOTYPE}.tmp
        """
        
##-------------------------------
# BOWTIE2 RIBO-RNA ALIGNMENT
##-------------------------------
rule rRNA_alignment:
    input:
        r1_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.repaired.fastq.gz',
        r2_repaired = OUTPUT+'{project}/sample_output/{sample}/{sample}_R2.repaired.fastq.gz'
    output:
        rRNA_bam = temp(OUTPUT+"{project}/sample_output/{sample}/{sample}_rRNA.bam"),
        rRNA_bwt_report = OUTPUT+"{project}/sample_output/{sample}/{sample}_bt2.output"
    params:
        rRNA_REF = MPATH + config["rRNA_"+GENOME_REF],
        BOWTIE2 = MPATH + config["BOWTIE2"]
    threads: ALIGN_THREADS
    log: "logs/{project}/{sample}.rRNA.log"
    shell:
        """
        {params.BOWTIE2} -x {params.rRNA_REF} \
            -1 {input.r1_repaired} -2 {input.r2_repaired} \
            -p {threads} -q \
            --dovetail \
            2> {output.rRNA_bwt_report} | samtools view -f3 -bh - > {output.rRNA_bam}
        """
        

##-------------------------------
# COMBINE OUTPUTS
##-------------------------------
rule output_stats:
    input:
        original_r1 = lambda w: get_fastq_path(w.sample) + '{sample}_R1.fastq.gz',
        clean_r1 = OUTPUT+'{project}/sample_output/{sample}/{sample}_R1.clean.fastq.gz',
        STAR_log = OUTPUT+'{project}/sample_output/{sample}/{sample}_Log.final.out',
        align_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
        dedup_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.DEDUPE.bam',
        ftcount_log = OUTPUT+'{project}/sample_output/{sample}/{sample}_feature_counts.txt.summary',
#         ftcount_log_NODEDUPE = OUTPUT+'{project}/sample_output/{sample}/{sample}_feature_counts.NONDEDUPE.txt.summary',
        exon_intron_stats = OUTPUT+'{project}/sample_output/{sample}/{sample}_exon_intron.stats',
        rRNA_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.RIBO.bam',
        MT_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.MT.bam',
        nucleosomal_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.NUCLEAR.bam',
        rnaseq_qc_results = OUTPUT + "{project}/sample_output/{sample}/rnaseq_qc_results.txt"
    output:
        stats = OUTPUT+'{project}/sample_output/{sample}/{sample}_mapping_stats.tsv'
    params:
        lib_prep = lambda w: sample_information.loc[w.sample, "prep_type"]
    shell:
        """
        echo -e sample_id'\t'{wildcards.sample} > {output.stats}
        
        original_reads=$(pigz -dc {input.original_r1} | wc -l)
        original_reads=$((original_reads / 4))
        echo -e reads_all'\t'${{original_reads}} >> {output.stats}
        
        
        trimmed_reads=$(pigz -dc {input.clean_r1} | wc -l)
        trimmed_reads=$((trimmed_reads / 4))
        echo -e reads_trimmed_filtered'\t'${{trimmed_reads}} >> {output.stats}
        
        
        align_input=$(grep "Number of input reads"  {input.STAR_log} | grep -Eo '[0-9]*')
        align_unique=$(grep "Uniquely mapped reads number"  {input.STAR_log} | grep -Eo '[0-9]*')
        align_multi=$(grep "Number of reads mapped to multiple loci"  {input.STAR_log} | grep -Eo '[0-9]*')
        align_tooshort=$(grep "Number of reads unmapped: too short"  {input.STAR_log} | grep -Eo '[0-9]*')
        echo -e align_input'\t'${{align_input}} >> {output.stats}
        echo -e align_unique'\t'${{align_unique}} >> {output.stats}
        echo -e align_multi'\t'${{align_multi}} >> {output.stats}
        echo -e align_tooshort'\t'${{align_tooshort}} >> {output.stats}


        align_out=$(samtools view -c {input.align_bam})
        dedup_out=$(samtools view -c {input.dedup_bam})
        duplication_rate=$(python -c "print(1 - ($dedup_out / $align_out))")
        echo -e duplication_rate'\t'${{duplication_rate}} >> {output.stats}
        
        
        ftcount=$(grep "Assigned" {input.ftcount_log} | grep -Eo "[0-9]*")
        echo -e ftcount'\t'${{ftcount}} >> {output.stats}
        
        
        intron_exon_ratio=$(python -c "print($intronic_reads / $exonic_reads)")
        exon_ratio=$(python -c "print($exonic_reads / ($intronic_reads + $exonic_reads + $intergenic_reads))")
        echo -e intron_exon_ratio'\t'${{intron_exon_ratio}} >> {output.stats}
        echo -e exon_ratio'\t'${{exon_ratio}} >> {output.stats}
        
        
        bias53=$(grep "5'-3' bias"  {input.rnaseq_qc_results} | cut -d "=" -f2 | xargs)
        echo -e bias53'\t'${{bias53}} >> {output.stats}
        
        
        rRNA_reads=$(samtools view -c {input.rRNA_bam})
        rRNA_reads=$((rRNA_reads/2))
        mtRNA_reads=$(samtools view -c {input.MT_bam})
        mtRNA_reads=$((mtRNA_reads/2))
        nuclear_reads=$(samtools view -c {input.nucleosomal_bam})
        nuclear_reads=$((nuclear_reads/2))
        echo -e rRNA_reads'\t'${{rRNA_reads}} >> {output.stats}
        echo -e mtRNA_reads'\t'${{mtRNA_reads}} >> {output.stats}
        echo -e nuclear_reads'\t'${{nuclear_reads}} >> {output.stats}
        
        
        exonic_reads=$(sed -n 1p {input.exon_intron_stats})
        exonic_reads=$(python -c "print($exonic_reads / 2)")
        intronic_reads=$(sed -n 2p {input.exon_intron_stats})
        intronic_reads=$(python -c "print($intronic_reads / 2)")
        intergenic_reads=$(sed -n 3p {input.exon_intron_stats})
        intergenic_reads=$(python -c "print($intergenic_reads / 2)")
        echo -e exonic_reads'\t'${{exonic_reads}} >> {output.stats}
        echo -e intronic_reads'\t'${{intronic_reads}} >> {output.stats}
        echo -e intergenic_reads'\t'${{intergenic_reads}} >> {output.stats}
        
        
        datamash transpose < {output.stats} > {output.stats}_T
        mv {output.stats}_T {output.stats}                
        """

#         ftcount_NODEDUPE=$(grep "Assigned" {input.ftcount_log_NODEDUPE} | grep -Eo "[0-9]*")
#        echo -e ftcount_NODEDUPE'\t'${{ftcount_NODEDUPE}} >> {output.stats}



####################
# ARCHIVED CODE !!!!
####################


# ##-------------------------------
# # LENGTH PROFILES
# # - This rule wil extract the length profiles for each bam file and plot them using a custom script Alex Cheng wrote awhile back.

# rule length_profiles:
#     input:
#         bam =  OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.NUCLEAR.bam'
#     output:
#         lengths=OUTPUT+'{project}/sample_output/{sample}/{sample}_sort_dedup.lengths.gz',
#         length_counts = OUTPUT+'{project}/sample_output/{sample}/{sample}_sort_dedup.lengths.counts',
#         pdf=OUTPUT+'{project}/sample_output/{sample}/{sample}_FragsHistogram.pdf'
#     conda:
#         "../envs/cfRNA.yaml"
#     threads: qc_threads
#     shell:
#         """
#         samtools view -F 0x100 {input.bam} | awk -v OFS='\t' '{{if ($9>0 && $9<1000) print $1,$9}}' | gzip -9 > {output.lengths}
#         zcat {output.lengths} | cut -f2 | sort | uniq -c > {output.length_counts}
#         ~/miniconda3/envs/cfrna/bin/Rscript workflow/scripts/length_profiles/HistogramFragmentLengths.R {output.length_counts} {output.pdf}
#         """

# rule intron_exon_bams:
#     input:
#         bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
#     output:
#         bam_bed = temp(OUTPUT+'{project}/sample_output/{sample}/{sample}_BAM.bed'),
#         stats_file = OUTPUT+'{project}/sample_output/{sample}/{sample}_exon_intron.stats',
#         exon_bed_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_exon.bed',
#         intron_bed_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_intron.bed',
#         intergenic_bed_bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_intergenic.bed'
#     params:
#         exon_bed= lambda w: config["EXONS_"+config['REF_DICT'][w.sample] ],
#         intron_bed= lambda w: config["INTRONS_"+config['REF_DICT'][w.sample] ],
#         intergenic_bed= lambda w: config["INTERGENIC_"+config['REF_DICT'][w.sample] ]
#     threads: 5
#     shell:
#         """
#         samtools view -bf 0x2 {input.bam} | samtools sort -n - | bedtools bamtobed -bedpe -mate1 -i stdin | awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $1, $6, $2, $7, $8, $10 }}' - | awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{if ($3 > $2) {{ print }}}}' - > {output.bam_bed}.TMP
        
#         samtools view -bf 0x2 {input.bam} | samtools sort -n - | bedtools bamtobed -bedpe -mate1 -i stdin | awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $1, $6, $2, $7, $8, $10 }}' - | awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{if ($3 < $2) {{ print }}}}' - | awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $1, $3, $2, $4, $5, $6 }}' - | cat - {output.bam_bed}.TMP > {output.bam_bed}
         
#         rm {output.bam_bed}.TMP
         
        
#         bedtools intersect -wa -s -c -f 1 -a {output.bam_bed} -b {params.exon_bed} > {output.exon_bed_bam}
        
#         bedtools intersect -v -s -a {output.bam_bed} -b {params.exon_bed} |  bedtools intersect -wa -s -c -f 1 -a stdin -b {params.intron_bed} > {output.intron_bed_bam}
        
#         bedtools intersect -v -s -a {output.bam_bed} -b {params.intron_bed} | bedtools intersect -v -s -a stdin -b {params.exon_bed} |  bedtools intersect -wa -s -c -f 1 -a stdin -b {params.intergenic_bed} > {output.intergenic_bed_bam}
          
#         sed -n '$=' {output.exon_bed_bam} > {output.stats_file}
#         sed -n '$=' {output.intron_bed_bam} >> {output.stats_file}      
#         sed -n '$=' {output.intergenic_bed_bam} >> {output.stats_file}      
#         """   
        
    
# rule samtools_stats:
#     input:
#         bam = OUTPUT + '{project}/sample_output/{sample}/{sample}_sort_dedup.bam',
#         bai = OUTPUT + '{project}/sample_output/{sample}/{sample}_sort_dedup.bam.bai'
#     output:
#         stats = OUTPUT + '{project}/sample_output/{sample}/{sample}_samtools-stats.txt',
#         idxstats = OUTPUT + '{project}/sample_output/{sample}/{sample}_samtools-idxstat.txt'
#     threads: qc_threads
#     shell:
#         """
#         samtools stats -@ {threads} {input.bam} > {output.stats}
#         samtools idxstats {input.bam} > {output.idxstats}
#         """

# rule samtools_flagstats:
#     input:
#         bam = OUTPUT + '{project}/sample_output/{sample}/{sample}_sort_dedup.bam'
#     output:
#         flagstats = OUTPUT + '{project}/sample_output/{sample}/{sample}_sort_dedup.flagstats.txt'
#     threads: qc_threads
#     shell:
#         """
#         samtools flagstat -@ {threads} {input.bam} >{output.flagstats}
#         """

# ##-------------------------------
# # Qualimap - basic summary stats

# rule qualimap:
#     input:
#         bam = OUTPUT+'{project}/sample_output/{sample}/{sample}_Aligned.sortedByCoord.out.bam',
#     output:
#         output_temp = OUTPUT + "{project}/sample_output/{sample}/rnaseq_qc_results.txt"
#     params:
#         outdir = OUTPUT+'{project}/sample_output/{sample}/',
#         ANNOTATION = lambda w: config["ANNOTATION_"+config['REF_DICT'][w.sample] ]
#     threads: 5
#     shell:
#         """
#         /programs/qualimap_v2.2.1/qualimap rnaseq \
#             -bam {input} \
#             -gtf {params.ANNOTATION} \
#             --outdir {params.outdir} \
#             --outfile {wildcards.sample}_qualimap \
#             --sequencing-protocol strand-specific-reverse \
#             --paired \
#             --java-mem-size=10G
#         """