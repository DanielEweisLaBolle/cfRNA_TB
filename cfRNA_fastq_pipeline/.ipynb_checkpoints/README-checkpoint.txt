##----------------------------##
#        cfRNA PIPELINE        #
##----------------------------##
VERSION: 2 - The Return of the Snake
DATE: 7/8/2022
CREATOR: Conor Loy

            / . .\
            \  ---<
             \  /
   __________/ /
-=:___________/



##----------------------------##
## Protocol to run the pipeline:
	
	 https://docs.google.com/document/d/1JLV8GDG0ilhHuiFf95ZbDU45JzXg-855-xdhdSndfto/edit?usp=sharing

Further information can be found below.
##----------------------------##


##----------------------------##
## Rules:
    1. DO NOT EDIT ANY FILE IN THIS DIRECTORY - if you want to modify the pipeline please copy to your workdir and make edits there (see below). The pipeline will execute from anywhere on cbsuvlaminck4!
    2. DO NOT DELETE ANY OUTPUT FILES - if there are issues please contact the CREATOR (see above)
    3. DO NOT DELETE OR MOVE ANYTHING IN THE "./prep_tables/" FOLDER - Please just copy from the prep_tables folder to "sequencing_prep.tsv" if you want to run something
    
    
##----------------------------##
## Overview:
    Plug and play cfRNA pipeline built with Snakemake. To process your samples create a sequencing prep table with your sample information and run a basic Snakemake command (eg. "Snakemake --cores 45 -k"). The outputs will be in the "output" folder underneath the project title you provided (see below for more information about file structure).
   
   
##----------------------------##
## Copying to your directory:
    To copy this pipeline to your directory run the COPY_PIPELINE.sh bash script using the following format: "./COPY_PIPELINE.sh <destination_folder>". This will copy the "takara_human_V2" directory to your destination, exluding outputs. YOU DO NOT NEED TO COPY ANY REFERENCES OR CHANGE PATHS IN THE CONFIG FILE TO THEN EXECUTE THE PIPELINE.
    Feel free to manually copy outputs to your workdir, but I thought it would be best to not have everyone have a copy of every output if it is not being used.
    
    
##----------------------------##    
## Sequencing Prep Table
    Create a sequencing prep table using the columns outlined below and save a copy to the "prep_tables/" folder with the project appended to the file name
    Columns:
     sample_id -> Sample id which will be used for naming; your fastqs should be named <sample_id>_R[1,2].fastq.gz
    project_id -> Name of the project for which the sample is a part of
     prep_type -> The library prep used for this sample (SMARTer_pico_v2 / SMARTer_pico_v3)
          path -> Absolute path to the location of the fastqs
      
      
##----------------------------##    
## File structure
    - All outputs can be found under the "output/" folder, organized by project
    - Individual sample outputs will be found here: "output/<project>/sample_output/<sample_id>/"
    - Within a "output/<project>" directory there are summary tables, including:
    
                       <project>_biotypes.tsv -> biotype counts of reads assigned to genes
    <project>_BP_celltypes.protein_coding.tsv -> cell type deconvolution estimates using protein coding genes
                   <project>_BP_celltypes.tsv -> cell type deconvolution estimates
                 <project>_feature_counts.tsv -> count matrix
                  <project>_mapping_stats.tsv -> QC and mapping statistics from the pipeline (intron/exon ratio, alignment rate, etc.)
    


##----------------------------##    
## Changes
    - Added the deconvolution to the pipeline (updated version of BayesPrism as well)
    - Removed the execution bash script "XcfRNA", and have moved everything within a Snakemake Command
    - Aggregation now occurs via snakemake rules
    - Minor changes to the *_mapping_stats.tsv file structure
    - Cleaned up syntax
    - Moved all paths to the config file
