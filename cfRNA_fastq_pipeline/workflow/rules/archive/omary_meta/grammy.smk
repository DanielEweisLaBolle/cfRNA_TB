rule grammy: # MAJOR MOD. DO NOT USE R2 (doubles the amount of reads, but the paired-end awareness makes R2 obsolete)
	input:
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa',
		tblat1 = 'sample_output/grammy/{sample}/{sample}.tblat.1'
	output:
		nonhumanfa_gz = temp('sample_output/grammy/{sample}/{sample}.fa.gz'),
		nonhumanfasta_gz = temp('sample_output/grammy/{sample}/{sample}.fasta.gz'),
		rdt = 'sample_output/grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/grammy/{sample}/{sample}.lld',
		btp = 'sample_output/grammy/{sample}/{sample}.btp',
		est = 'sample_output/grammy/{sample}/{sample}.est',
		gra = 'sample_output/grammy/{sample}/{sample}.gra',
		avl = 'sample_output/grammy/{sample}/{sample}.avl'
	resources: mem_mb=1
	shell: #used to be cat R1 R2 | sed ....
		"""
		cat {input.nonhumanfa_r1} | sed 's/_1.*/-1/g' |  sed 's/_2.*/-2/g' | gzip -1 > {output.nonhumanfa_gz}
		cd sample_output/grammy/{wildcards.sample}
		python2.7 {GRAMMY_RDT} -t illumina . .
		python2.7 {GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {GRAMMY_REF_FASTA}
		python2.7 {GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
		python2.7 {GRAMMY_POST} {wildcards.sample}.est {GRAMMY_REF_FASTA} {wildcards.sample}.btp
		cd ../../../
		"""

rule annotate_grammy:
	input:
		rdt = 'sample_output/grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/grammy/{sample}/{sample}.lld',
		btp = 'sample_output/grammy/{sample}/{sample}.btp',
		est = 'sample_output/grammy/{sample}/{sample}.est',
		gra = 'sample_output/grammy/{sample}/{sample}.gra',
		avl = 'sample_output/grammy/{sample}/{sample}.avl',
		stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		tblat1 = 'sample_output/grammy/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/grammy/{sample}/{sample}.tab',
		anno = 'sample_output/grammy/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='grammy/{sample}/',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
		Rscript scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
		Rscript scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {input.stats} {input.LUT} {output.anno}
		"""

rule grammy_CFS: # For CFS reads, the identified contaminants are removed from the nonhuman FASTAs
	input:
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa',
		tblat1 = 'sample_output/filtered_grammy/{sample}/{sample}.tblat.1',
		removed_reads = 'sample_output/filtered_grammy/{sample}/removed.{sample}.blast'
	output:
		to_remove = temp('sample_output/filtered_grammy/{sample}/{sample}.toremove'),
		to_keep = temp('sample_output/filtered_grammy/{sample}/{sample}.tokeep'),
		to_remove2 = temp('sample_output/filtered_grammy/{sample}/{sample}.toremove2'),
		nonhumanfa_gz = temp('sample_output/filtered_grammy/{sample}/{sample}.fa.gz'),
		nonhumanfasta_gz = temp('sample_output/filtered_grammy/{sample}/{sample}.fasta.gz'),
		rdt = 'sample_output/filtered_grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/filtered_grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/filtered_grammy/{sample}/{sample}.lld',
		btp = 'sample_output/filtered_grammy/{sample}/{sample}.btp',
		est = 'sample_output/filtered_grammy/{sample}/{sample}.est',
		gra = 'sample_output/filtered_grammy/{sample}/{sample}.gra',
		avl = 'sample_output/filtered_grammy/{sample}/{sample}.avl'
	resources: mem_mb=1
	shell: #used to be cat R1 R2 | sed ....
		"""
		tail -n+1 {input.removed_reads} | cut -f1 | sort -u > {output.to_remove}
		cut -f1 {input.tblat1} | sort -u > {output.to_keep}
		echo "1"
		comm -23 {output.to_remove} {output.to_keep} > {output.to_remove2}
		echo "2"
		cat {input.nonhumanfa_r1} | sed 's/_1.*/-1/g' |  sed 's/_2.*/-2/g' | \
			./software/bbmap/filterbyname.sh in=stdin.fasta out={output.nonhumanfa_gz} names={output.to_remove2} fastawrap=500
		cd sample_output/filtered_grammy/{wildcards.sample}
		python2.7 {GRAMMY_RDT} -t illumina . .
		python2.7 {GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {GRAMMY_REF_FASTA}
		python2.7 {GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
		python2.7 {GRAMMY_POST} {wildcards.sample}.est {GRAMMY_REF_FASTA} {wildcards.sample}.btp
		cd ../../../
		"""

rule annotate_grammy_CFS:
	input:
		rdt = 'sample_output/filtered_grammy/{sample}/{sample}.rdt',
		mtx = 'sample_output/filtered_grammy/{sample}/{sample}.mtx',
		lld = 'sample_output/filtered_grammy/{sample}/{sample}.lld',
		btp = 'sample_output/filtered_grammy/{sample}/{sample}.btp',
		est = 'sample_output/filtered_grammy/{sample}/{sample}.est',
		gra = 'sample_output/filtered_grammy/{sample}/{sample}.gra',
		avl = 'sample_output/filtered_grammy/{sample}/{sample}.avl',
		stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		tblat1 = 'sample_output/filtered_grammy/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/filtered_grammy/{sample}/{sample}.tab',
		anno = 'sample_output/filtered_grammy/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='filtered_grammy/{sample}/',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
		Rscript scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
		Rscript scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {input.stats} {input.LUT} {output.anno}
		"""

rule grammy_CFS2: # For CFS reads, the identified contaminants are removed from the nonhuman FASTAs
	input:
		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa',
		tblat1 = 'sample_output/filtered_grammy/{sample}/{sample}.tblat.1',
		removed_reads = 'sample_output/filtered_grammy/{sample}/removed.{sample}.blast'
	output:
		tblat1 = 'sample_output/filtered_grammy2/{sample}/{sample}.tblat.1',
		removed_reads = 'sample_output/filtered_grammy2/{sample}/removed.{sample}.blast',
		to_keep = temp('sample_output/filtered_grammy2/{sample}/{sample}.toremove'),
		nonhumanfa_gz = temp('sample_output/filtered_grammy2/{sample}/{sample}.fa.gz'),
		nonhumanfasta_gz = temp('sample_output/filtered_grammy2/{sample}/{sample}.fasta.gz'),
		rdt = 'sample_output/filtered_grammy2/{sample}/{sample}.rdt',
		mtx = 'sample_output/filtered_grammy2/{sample}/{sample}.mtx',
		lld = 'sample_output/filtered_grammy2/{sample}/{sample}.lld',
		btp = 'sample_output/filtered_grammy2/{sample}/{sample}.btp',
		est = 'sample_output/filtered_grammy2/{sample}/{sample}.est',
		gra = 'sample_output/filtered_grammy2/{sample}/{sample}.gra',
		avl = 'sample_output/filtered_grammy2/{sample}/{sample}.avl'
	resources: mem_mb=1
	shell: #used to be cat R1 R2 | sed ....
		"""
		cp {input.tblat1} {output.tblat1}
		cp {input.removed_reads} {output.removed_reads}
		cut -f1 {output.tblat1} | sort -u > {output.to_keep}
		cat {input.nonhumanfa_r1} | sed 's/_1.*/-1/g' |  sed 's/_2.*/-2/g' | \
			./software/bbmap/filterbyname.sh include=t in=stdin.fasta out={output.nonhumanfa_gz} names={output.to_keep} fastawrap=500
		cd sample_output/filtered_grammy2/{wildcards.sample}
		python2.7 {GRAMMY_RDT} -t illumina . .
		python2.7 {GRAMMY_PRE} -q "40,75,-5" {wildcards.sample} {GRAMMY_REF_FASTA}
		python2.7 {GRAMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx
		python2.7 {GRAMMY_POST} {wildcards.sample}.est {GRAMMY_REF_FASTA} {wildcards.sample}.btp
		cd ../../../
		"""

rule annotate_grammy_CFS2:
	input:
		rdt = 'sample_output/filtered_grammy2/{sample}/{sample}.rdt',
		mtx = 'sample_output/filtered_grammy2/{sample}/{sample}.mtx',
		lld = 'sample_output/filtered_grammy2/{sample}/{sample}.lld',
		btp = 'sample_output/filtered_grammy2/{sample}/{sample}.btp',
		est = 'sample_output/filtered_grammy2/{sample}/{sample}.est',
		gra = 'sample_output/filtered_grammy2/{sample}/{sample}.gra',
		avl = 'sample_output/filtered_grammy2/{sample}/{sample}.avl',
		stats = 'sample_output/stats/{sample}_mapping_stats.txt',
		tblat1 = 'sample_output/filtered_grammy2/{sample}/{sample}.tblat.1',
		LUT = 'LUTGrammy/taxids_names_lengths_tax.tab'
	output:
		tab='sample_output/filtered_grammy2/{sample}/{sample}.tab',
		anno = 'sample_output/filtered_grammy2/{sample}/{sample}.grammy.tab'
	params:
		DIR='sample_output/',
		DB='filtered_grammy2/{sample}/',
		LUT='LUTGrammy/taxids_names_lengths_tax.tab'
	shell:
		"""
		Rscript scripts/metagenome/filter_gra_file.R {input.gra} {output.tab} {wildcards.sample}
		Rscript scripts/metagenome/annotate_grammy_apc.R {output.tab} {input.tblat1} {input.stats} {input.LUT} {output.anno}
		"""

# rule cleanup:
# 	input:
# 		anno = 'sample_output/grammy/{sample}/{sample}.grammy.tab',
# 		decon_r1 = 'sample_output/decontaminate/{sample}_decontaminated_R1.fastq',
# 		decon_r2 = 'sample_output/decontaminate/{sample}_decontaminated_R2.fastq',
# 		nonhumanfa_r1 = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa',
# 		nonhumanfa_r2 = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa'
# 	output:
# 		decon_r1_gz = 'sample_output/decontaminate/{sample}_decontaminated_R1.fastq.gz',
# 		decon_r2_gz = 'sample_output/decontaminate/{sample}_decontaminated_R2.fastq.gz',
# 		nonhumanfa_r1_gz = 'sample_output/nonhuman_fasta/R1/{sample}_R1.fa.gz',
# 		nonhumanfa_r2_gz = 'sample_output/nonhuman_fasta/R2/{sample}_R2.fa.gz'
# 	shell:
# 		"""
# 		gzip -c {input.decon_r1} > {output.decon_r1_gz}
# 		gzip -c {input.decon_r2} > {output.decon_r2_gz}
# 		gzip -c {input.nonhumanfa_r1} > {output.nonhumanfa_r1_gz}
# 		gzip -c {input.nonhumanfa_r2} > {output.nonhumanfa_r2_gz}
# 		"""
