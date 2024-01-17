rm(list=ls())

##------------------------------------------

library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)

'%ni%' <- Negate('%in%')

##------------------------------------------

args = commandArgs(trailingOnly=TRUE)
project <- args[1]


output_dir = paste0("./output/",project)

seq_df = read.delim("./prep_tables/sequencing_prep_ALL.tsv")
all_IDs <- seq_df[seq_df$project_id == project,'sample_id']


##------------------------------------------
# rRNA CONTAMINATION

rRNA <- function(ID,output_dir){
return(as.numeric(gsub("\\%","",unlist(strsplit(unlist(strsplit(readLines(paste0(output_dir,"/",ID,"/",ID,"_bt2.output"))[4],"\\("))[2],"\\)"))[1])))
}


##------------------------------------------
# DNA CONTAMINATION

DNA <- function(ID, output_dir){
df = read.delim(paste0(output_dir,"/",ID,"/",ID,".metrics.tsv"))
return(as.numeric(df[df$Sample == 'Intronic Rate',2]) / as.numeric(df[df$Sample == 'Exonic Rate',2]))
}

##------------------------------------------
# RNA DEGREDATION
# RNA_DEG <- function(){
	
# }

##------------------------------------------
##------------------------------------------
# COMBINE ALL

rRNA_rates = unlist(lapply(all_IDs,rRNA,output_dir))
DNA_rates = unlist(lapply(all_IDs,DNA,output_dir))

summary_df = data.frame('sample_id' = all_IDs,
						'rRNA_rate' = rRNA_rates,
					   'intron_exon_ratio' = DNA_rates 
					   )

write.table(summary_df,paste0(output_dir,"/",project,"_summary-stats.tsv"),
		   sep = "\t", quote = FALSE,row.names=FALSE)


# TO ADD:
# - biotype distribution
# intron / exon over lap plots