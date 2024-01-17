

rule hs_blastn_stds:
	input:
		nonhumanfa = 'sample_output/nonhuman_fasta_stds/{read}/{sample}_{read}.fa',
		db = '/fs/cbsuvlaminck/workdir/omm35/WGBS_pipeline/databases/blast/stds/NCBIGenomes06.fna',
		gi_to_taxid = '/fs/cbsuvlaminck/workdir/omm35/WGBS_pipeline/databases/blast/NCBIGenomes06.gis.taxids',
		grammy_db_proof = '/fs/cbsuvlaminck/workdir/omm35/WGBS_pipeline/databases/logs/grammy/stds_grammy_prompts',
		obinary = '/fs/cbsuvlaminck/workdir/omm35/WGBS_pipeline/databases/blast/stds/NCBIGenomes06.fna.counts.obinary'
	output:
		blast_outfmt6 = temp('sample_output/blast/stds/{read}/{sample}.{read}.outfmt6'),
		human_like = temp('sample_output/blast/stds/{read}/{sample}.{read}.humanlike'),
		rejected = temp('sample_output/blast/stds/{read}/{sample}.{read}.rejected')
	threads: blast_threads
	resources:
		mem_mb=20000
	shell:
		"""
		{HSBLASTN} align -query {input.nonhumanfa} \
                    	-db {input.db} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 | \
		python scripts/metagenome/get_taxid_filter_strand.py --filename_out {output.blast_outfmt6} \
															--gi_to_tax {input.gi_to_taxid} \
															--conversion stds \
															--human_like {output.human_like} \
															--rejected_hits {output.rejected}
		"""
								#-window_masker_db {input.obinary} \

rule compress_stds_outfmt6:
	input:
		outfmt6 = 'sample_output/blast/stds/{read}/{sample}.{read}.outfmt6'
	output:
		gzip = 'sample_output/blast/stds/{read}/{sample}.{read}.outfmt6.gz'
	threads: 10
	shell:
		"""
		pigz --best -c -p {threads} {input.outfmt6} > {output.gzip}
		"""

rule filter_blastn_stds:
	input:
		R1 = 'sample_output/blast/stds/R1/{sample}.R1.outfmt6.gz',
		R2 = 'sample_output/blast/stds/R2/{sample}.R2.outfmt6.gz',
		taxid_lengths = '/fs/cbsuvlaminck/workdir/omm35/WGBS_pipeline/databases/GenomeDB/taxids_lengths.txt'
	output:
		tblatpe= temp('sample_output/blast/consolidated_blast_stds/{sample}.tblat.pe'),
		tblat1 = temp('sample_output/grammy_stds/{sample}/{sample}.tblat.1')
	threads:
		10
	shell:
		"""
		filesize=$(du {input.R1} | cut -f1)
		if [ $filesize -gt 3696400 ]
		then
			bash scripts/metagenome/testing_split2.sh {input.R1} {input.R2} {wildcards.sample}_STDS {input.taxid_lengths} {threads} {output.tblatpe}

		else
			Rscript scripts/metagenome/filter_paired_end_blast2.R {input.R1} {input.R2} {output.tblatpe} {input.taxid_lengths}
		fi
		cat {output.tblatpe} | grep -v qseqid | awk '{{print $2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$3,$1}}' | tr ' ' '\t' > {output.tblat1}
		"""


def aggregate_tblat_files(wildcards):
	"""
	we have this function to figure out which inputs are hg19, bisulfite, etc.
	"""
	sample_name, prep_type, seq_mode, abundance_control, sample_type, seq_type = get_sample_info(wildcards)

	if "2x" in seq_mode and seq_type == "bisulfite":
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta_wgbs/R1/{sample}_R1.fa'
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta_wgbs/R2/{sample}_R2.fa'
		tblatpe = 'sample_output/blast/consolidated_blast_wgbs/{sample}.tblat.pe'
		tblat1 = 'sample_output/grammy_wgbs/{sample}/{sample}.tblat.1'

	if "2x" in seq_mode and seq_type == "standard":
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta_stds/R1/{sample}_R1.fa'
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta_stds/R2/{sample}_R2.fa'
		tblatpe = 'sample_output/blast/consolidated_blast_stds/{sample}.tblat.pe'
		tblat1 = 'sample_output/grammy_stds/{sample}/{sample}.tblat.1'

	return[nonhumanfa_r1, nonhumanfa_r2, tblatpe, tblat1]

rule aggregate_tblats:
	input:
		aggregate_tblat_files
	output:
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa',
		tblatpe = 'sample_output/blast/consolidated_blast/{sample}.tblat.pe',
		tblat1 = 'sample_output/grammy/{sample}/{sample}.tblat.1'
	shell:
		"""
		cp {input[0]} {output.nonhumanfa_r1}
		cp {input[1]} {output.nonhumanfa_r2}
		cp {input[2]} {output.tblatpe}
		cp {input[3]} {output.tblat1}
		"""
