
rule trim_unmapped:
    input: 
        r1 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R1.fastq',
        r2 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R2.fastq'
    output:
        r1_trim = OUTPUT+'{project}/{sample}/{sample}_unmapped_R1_trimmed.fq',
        r2_trim = OUTPUT+'{project}/{sample}/{sample}_unmapped_R2_trimmed.fq',
        r1_unpaired = OUTPUT+'{project}/{sample}/{sample}_unmapped_R1_unpaired_trimmed.fq',
        r2_unpaired = OUTPUT+'{project}/{sample}/{sample}_unmapped_R2_unpaired_trimmed.fq',
        trim_log = OUTPUT+'{project}/{sample}/{sample}_unmapped_trimmomcatic.log'
    threads: trim_threads
    params:
        output_dir = OUTPUT+"{project}/{sample}/",
        adapters = config["TRUSEQ_ADAPT"]
    log: 'logs/{project}/{sample}.trim.log'
    shell:
        """
        java -jar /programs/trimmomatic/trimmomatic-0.39.jar \
            PE -phred33 -threads {threads} -trimlog {output.trim_log}\
            {input.r1} {input.r2} \
            {output.r1_trim} {output.r1_unpaired} \
            {output.r2_trim} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:36
        """      
        
        
rule decontaminate_stds:
    input:
        unmapped_R1 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R1_trimmed.fq',
        unmapped_R2 = OUTPUT+'{project}/{sample}/{sample}_unmapped_R2_trimmed.fq',
        both = config["BOTH_REF"] 
    output:
        r1_pass1 = temp(OUTPUT+'{project}/{sample}/{sample}.pass1_R1.fastq'),
        r2_pass1 = temp(OUTPUT+'{project}/{sample}/{sample}.pass1_R2.fastq'),
        decon_r1 = temp(OUTPUT+'{project}/{sample}/{sample}_decontaminated_R1.fastq'),
        decon_r2 = temp(OUTPUT+'{project}/{sample}/{sample}_decontaminated_R2.fastq'),
        nonhumanfa_r1 = temp(OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1.fa'),
        nonhumanfa_r2 = temp(OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2.fa')
    resources:
        mem_mb=80000
    params:
        BBDUK = config["BBDUK"],
        SEQTK = config["SEQTK"],
        genome = lambda w: config["REF_GENOME_"+config['REF_DICT'][w.sample]] 
    shell:
        """
        {params.BBDUK} in1={input.unmapped_R1} in2={input.unmapped_R2} \
                out1={output.r1_pass1} out2={output.r2_pass1} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={params.genome} k=50
        {params.BBDUK} in1={output.r1_pass1} in2={output.r2_pass1} \
                out1={output.decon_r1} out2={output.decon_r2} \
                -Xmx{resources.mem_mb}m \
                prealloc=t \
                ref={input.both} k=20
        fastq_to_fasta -Q33 -i {output.decon_r1} -o {output.nonhumanfa_r1}
        fastq_to_fasta -Q33 -i {output.decon_r2} | {params.SEQTK} seq -r - > {output.nonhumanfa_r2}
        """

        
###################################################################################
#BLAST
###################################################################################
rule hs_blastn_stds_R1:
    input:
        nonhumanfa = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1.fa',
        db = '/workdir/omm35/WGBS_pipeline/databases/blast/stds/NCBIGenomes06.fna',
        gi_to_taxid = '/workdir/omm35/WGBS_pipeline/databases/blast/NCBIGenomes06.gis.taxids',
        grammy_db_proof = '/workdir/omm35/WGBS_pipeline/databases/logs/grammy/stds_grammy_prompts',
        obinary = '/workdir/omm35/WGBS_pipeline/databases/blast/stds/NCBIGenomes06.fna.counts.obinary'
    output:
        blast_outfmt6 = temp(OUTPUT+'{project}/{sample}/blast/{sample}.R1.outfmt6'),
        human_like = temp(OUTPUT+'{project}/{sample}/{sample}.R1.humanlike'),
        rejected = temp(OUTPUT+'{project}/{sample}/blast/{sample}.R1.rejected')
    threads: blast_threads
    resources:
        mem_mb=20000
    params:
        HSBLASTN = config["HSBLASTN"]
    shell:
        """
        {params.HSBLASTN} align -query {input.nonhumanfa} \
                        -db {input.db} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 | \
        python workflow/scripts/metagenomics_scripts/get_taxid_filter_strand.py --filename_out {output.blast_outfmt6} \
                                                            --gi_to_tax {input.gi_to_taxid} \
                                                            --conversion stds \
                                                            --human_like {output.human_like} \
                                                            --rejected_hits {output.rejected}
        """
        
rule compress_stds_outfmt6_R1:
    input:
        outfmt6 = OUTPUT+'{project}/{sample}/blast/{sample}.R1.outfmt6'
    output:
        gzip = OUTPUT+'{project}/{sample}/blast/{sample}.R1.outfmt6.gz'
    threads: 10
    shell:
        """
        pigz --best -c -p {threads} {input.outfmt6} > {output.gzip}
        """
        
rule hs_blastn_stds_R2:
    input:
        nonhumanfa = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2.fa',
        db = '/workdir/omm35/WGBS_pipeline/databases/blast/stds/NCBIGenomes06.fna',
        gi_to_taxid = '/workdir/omm35/WGBS_pipeline/databases/blast/NCBIGenomes06.gis.taxids',
        grammy_db_proof = '/workdir/omm35/WGBS_pipeline/databases/logs/grammy/stds_grammy_prompts',
        obinary = '/workdir/omm35/WGBS_pipeline/databases/blast/stds/NCBIGenomes06.fna.counts.obinary'
    output:
        blast_outfmt6 = temp(OUTPUT+'{project}/{sample}/blast/{sample}.R2.outfmt6'),
        human_like = temp(OUTPUT+'{project}/{sample}/blast/{sample}.R2.humanlike'),
        rejected = temp(OUTPUT+'{project}/{sample}/blast/{sample}.R2.rejected')
    threads: blast_threads
    params:
        HSBLASTN = config["HSBLASTN"]
    resources:
        mem_mb=20000
    shell:
        """
        {params.HSBLASTN} align -query {input.nonhumanfa} \
                        -db {input.db} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 | \
        python workflow/scripts/metagenomics_scripts/get_taxid_filter_strand.py --filename_out {output.blast_outfmt6} \
                                                            --gi_to_tax {input.gi_to_taxid} \
                                                            --conversion stds \
                                                            --human_like {output.human_like} \
                                                            --rejected_hits {output.rejected}
        """
                            
rule compress_stds_outfmt6_R2:
    input:
        outfmt6 = OUTPUT+'{project}/{sample}/blast/{sample}.R2.outfmt6'
    output:
        gzip = OUTPUT+'{project}/{sample}/blast/{sample}.R2.outfmt6.gz'
    threads: 10
    shell:
        """
        pigz --best -c -p {threads} {input.outfmt6} > {output.gzip}
        """
        


rule filter_blastn_stds:
    input:
        R1 = OUTPUT+'{project}/{sample}/blast/{sample}.R1.outfmt6.gz',
        R2 = OUTPUT+'{project}/{sample}/blast/{sample}.R2.outfmt6.gz',
        taxid_lengths = '/workdir/omm35/WGBS_pipeline/databases/GenomeDB/taxids_lengths.txt'
    output:
        tblatpe= temp(OUTPUT+'{project}/{sample}/grammy/{sample}.tblat.pe'),
        tblat1 = temp(OUTPUT+'{project}/{sample}/grammy/{sample}.tblat.1')
    threads:
        10
    shell:
        """
        filesize=$(du {input.R1} | cut -f1)
        if [ $filesize -gt 3696400 ]
        then
            bash workflow/scripts/metagenomics_scripts/testing_split2.sh {input.R1} {input.R2} {wildcards.sample}_STDS {input.taxid_lengths} {threads} {output.tblatpe}

        else
            Rscript workflow/scripts/metagenomics_scripts/filter_paired_end_blast2.R {input.R1} {input.R2} {output.tblatpe} {input.taxid_lengths}
        fi
        cat {output.tblatpe} | grep -v qseqid | awk '{{print $2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$3,$1}}' | tr ' ' '\t' > {output.tblat1}
        """

###################################################################################
#GRAMMY
###################################################################################
        
rule grammy: # MAJOR MOD. DO NOT USE R2 (doubles the amount of reads, but the paired-end awareness makes R2 obsolete)
    input:
        nonhumanfa_r1 = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R1.fa',
        nonhumanfa_r2 = OUTPUT+'{project}/{sample}/{sample}_nonhuman_R2.fa',
        tblat1 = OUTPUT+'{project}/{sample}/grammy/{sample}.tblat.1'
    output:
        nonhumanfa_gz = temp(OUTPUT+'{project}/{sample}/grammy/{sample}.fa.gz'),
        nonhumanfasta_gz = temp(OUTPUT+'{project}/{sample}/grammy/{sample}.fasta.gz'),
        rdt = OUTPUT+'{project}/{sample}/grammy/{sample}.rdt',
        mtx = OUTPUT+'{project}/{sample}/grammy/{sample}.mtx',
        lld = OUTPUT+'{project}/{sample}/grammy/{sample}.lld',
        btp = OUTPUT+'{project}/{sample}/grammy/{sample}.btp',
        est = OUTPUT+'{project}/{sample}/grammy/{sample}.est',
        gra = OUTPUT+'{project}/{sample}/grammy/{sample}.gra',
        avl = OUTPUT+'{project}/{sample}/grammy/{sample}.avl'
    resources: mem_mb=1
    params:
        output_dir = OUTPUT+'{project}/{sample}/grammy/', #output_dir = OUTPUT+'{project}/{sample}/grammy/',
        GRAMMY_RDT = config['GRAMMY_RDT'],
        GRAMMY_PRE = config['GRAMMY_PRE'],
        GRAMMY_EM = config['GRAMMY_EM'],
        GRAMMY_POST = config['GRAMMY_POST'],
        GRAMMY_REF_FASTA = config['GRAMMY_REF_FASTA']
    shell: #used to be cat R1 R2 | sed ....
        """
        cat {input.nonhumanfa_r1} | sed 's/_1.*/-1/g' |  sed 's/_2.*/-2/g' | gzip -1 > {output.nonhumanfa_gz}
        cd {params.output_dir}
        python2.7 {params.GRAMMY_RDT} -t illumina . .
        python2.7 {params.GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {params.GRAMMY_REF_FASTA}
        python2.7 {params.GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
        python2.7 {params.GRAMMY_POST} {wildcards.sample}.est {params.GRAMMY_REF_FASTA} {wildcards.sample}.btp
        cd ../../../
        """
# python2.7
rule annotate_grammy:
    input:
        rdt = OUTPUT+'{project}/{sample}/grammy/{sample}.rdt',
        mtx = OUTPUT+'{project}/{sample}/grammy/{sample}.mtx',
        lld = OUTPUT+'{project}/{sample}/grammy/{sample}.lld',
        btp = OUTPUT+'{project}/{sample}/grammy/{sample}.btp',
        est = OUTPUT+'{project}/{sample}/grammy/{sample}.est',
        gra = OUTPUT+'{project}/{sample}/grammy/{sample}.gra',
        avl = OUTPUT+'{project}/{sample}/grammy/{sample}.avl',
        stats = OUTPUT+'{project}/{sample}/{sample}_mapping_stats.txt',
        tblat1 = OUTPUT+'{project}/{sample}/grammy/{sample}.tblat.1',
        LUT = config['LUTGrammy']
    output:
        tab=OUTPUT+'{project}/{sample}/grammy/{sample}.tab',
        anno = OUTPUT+'{project}/{sample}/grammy/{sample}.grammy.tab'
#     params:
#         DIR=OUTPUT+'{project}/{sample}/grammy/,
#         DB='grammy/{sample}/',
#         LUT='LUTGrammy/taxids_names_lengths_tax.tab'
    shell:
        """
        Rscript workflow/scripts/metagenomics_scripts/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
        Rscript workflow/scripts/metagenomics_scripts/annotate_grammy_apc.R {output.tab} {input.tblat1} {input.stats} {input.LUT} {output.anno}
        """