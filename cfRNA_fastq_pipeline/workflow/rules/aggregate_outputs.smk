##-------------------------------
# COUNT MATRIX
##-------------------------------
def getSamplesInProject__agg_feature_counts(wildcards):
    return [ OUTPUT + wildcards.project + '/sample_output/'+s+'/'+s+'_FTCOUNTS.txt' for s in samplesInProjects[wildcards.project]]

def getSamplesInProject__agg_feature_counts_BIOTYPES(wildcards):
    return [ OUTPUT + wildcards.project + '/sample_output/'+s+'/'+s+'_BIOTYPES.txt' for s in samplesInProjects[wildcards.project]]

rule agg_feature_counts:
    input:
        ftcounts = getSamplesInProject__agg_feature_counts,
        biotypes = getSamplesInProject__agg_feature_counts_BIOTYPES
    output:
        agg_ftcounts= OUTPUT + '{project}/{project}_feature_counts.tsv',
        agg_biotypes = OUTPUT + '{project}/{project}_biotypes.tsv'
    threads: DECON_THREADS
    shell:
        """
        paste -d " " {input.ftcounts} > {output.agg_ftcounts}.tmp
        cut -d " " -f 1 {output.agg_ftcounts}.tmp > {output.agg_ftcounts}.tmp.geneids
        awk '{{for (i=2;i<=NF;i+=2) printf $i OFS}} {{print ""}}' {output.agg_ftcounts}.tmp > {output.agg_ftcounts}.tmp.cnts
        paste -d " " {output.agg_ftcounts}.tmp.geneids {output.agg_ftcounts}.tmp.cnts | awk -v OFS="\t" '$1=$1' - > {output.agg_ftcounts}
        rm {output.agg_ftcounts}.tmp*
        
        paste -d " " {input.biotypes} > {output.agg_biotypes}.tmp
        cut -d " " -f 1 {output.agg_biotypes}.tmp > {output.agg_biotypes}.tmp.biotypes
        awk '{{for (i=2;i<=NF;i+=2) printf $i OFS}} {{print ""}}' {output.agg_biotypes}.tmp > {output.agg_biotypes}.tmp.cnts
        paste -d " " {output.agg_biotypes}.tmp.biotypes {output.agg_biotypes}.tmp.cnts | awk -v OFS="\t" '$1=$1' - > {output.agg_biotypes}
        rm {output.agg_biotypes}.tmp*
        
        """
#         paste {input.ftcounts} > {output.agg_ftcounts}
#         awk '{for (i=2;i<=NF;i+=2) printf $i OFS} {print ""}' {output.agg_ftcounts}.tmp | {output.agg_ftcounts}
#         awk '{for(x=1;x<=NF;x++)if(x % 2)printf "%s", $x (x == NF || x == (NF-1)?"\n":" ")}' ahj.txt

##-------------------------------
# DECONVOLUTION
##-------------------------------       
def getSamplesInProject__agg_BP_decon(wildcards):
    return [ OUTPUT + wildcards.project + '/sample_output/'+s+'/'+s+'_BP_celltypes.tsv' for s in samplesInProjects[wildcards.project]]

def getSamplesInProject__agg_BP_decon_pc(wildcards):
    return [ OUTPUT + wildcards.project + '/sample_output/'+s+'/'+s+'_BP_celltypes.protein_coding.tsv' for s in samplesInProjects[wildcards.project]]

def getSamplesInProject__agg_BP_decon_pc_sig(wildcards):
    return [ OUTPUT + wildcards.project + '/sample_output/'+s+'/'+s+'_BP_celltypes.protein_coding.marker_genes.tsv' for s in samplesInProjects[wildcards.project]]

rule agg_BP_decon:
    input:
        celltype_frac = getSamplesInProject__agg_BP_decon,
        celltype_frac_pc = getSamplesInProject__agg_BP_decon_pc,
        celltype_frac_pc_sig = getSamplesInProject__agg_BP_decon_pc_sig,
    output:
        agg_celltype_frac = OUTPUT + '{project}/{project}_BP_celltypes.tsv',
        agg_celltype_frac_pc = OUTPUT + '{project}/{project}_BP_celltypes.protein_coding.tsv',
        agg_celltype_frac_pc_sig = OUTPUT + '{project}/{project}_BP_celltypes.protein_coding.marker_genes.tsv'
    threads: DECON_THREADS
    shell:
        """
        cat {input.celltype_frac} > {output.agg_celltype_frac}.tmp
        head -n 1 {output.agg_celltype_frac}.tmp > {output.agg_celltype_frac}.tmp.header
        sed -n 'n;p' {output.agg_celltype_frac}.tmp > {output.agg_celltype_frac}.tmp.body
        cat {output.agg_celltype_frac}.tmp.header {output.agg_celltype_frac}.tmp.body > {output.agg_celltype_frac}
        rm {output.agg_celltype_frac}.tmp*
        
        
        cat {input.celltype_frac_pc} > {output.agg_celltype_frac_pc}.tmp
        head -n 1 {output.agg_celltype_frac_pc}.tmp > {output.agg_celltype_frac_pc}.tmp.header
        sed -n 'n;p' {output.agg_celltype_frac_pc}.tmp > {output.agg_celltype_frac_pc}.tmp.body
        cat {output.agg_celltype_frac_pc}.tmp.header {output.agg_celltype_frac_pc}.tmp.body > {output.agg_celltype_frac_pc}
        rm {output.agg_celltype_frac_pc}.tmp*
        
        cat {input.celltype_frac_pc_sig} > {output.agg_celltype_frac_pc_sig}.tmp
        head -n 1 {output.agg_celltype_frac_pc_sig}.tmp > {output.agg_celltype_frac_pc_sig}.tmp.header
        sed -n 'n;p' {output.agg_celltype_frac_pc_sig}.tmp > {output.agg_celltype_frac_pc_sig}.tmp.body
        cat {output.agg_celltype_frac_pc_sig}.tmp.header {output.agg_celltype_frac_pc_sig}.tmp.body > {output.agg_celltype_frac_pc_sig}
        rm {output.agg_celltype_frac_pc_sig}.tmp*
        """
    
    
##-------------------------------
# MAPPING STATS
##-------------------------------       
def getSamplesInProject__agg_output_stats(wildcards):
    return [ OUTPUT + wildcards.project + '/sample_output/'+s+'/'+s+'_mapping_stats.tsv' for s in samplesInProjects[wildcards.project]]


rule agg_output_stats:
    input:
        getSamplesInProject__agg_output_stats
#         stats = expand( OUTPUT + '{{project}}/sample_output/{sample}/{sample}_mapping_stats.tsv', sample = SAMPLES )
    output:
        agg_stats = OUTPUT+'{project}/{project}_mapping_stats.tsv'
    shell:
        """
        cat {input} > {output.agg_stats}.tmp
        head -n 1 {output.agg_stats}.tmp > {output.agg_stats}.tmp.header
        sed -n 'n;p' {output.agg_stats}.tmp > {output.agg_stats}.tmp.body
        cat {output.agg_stats}.tmp.header {output.agg_stats}.tmp.body > {output.agg_stats}
        rm {output.agg_stats}.tmp*
        """
