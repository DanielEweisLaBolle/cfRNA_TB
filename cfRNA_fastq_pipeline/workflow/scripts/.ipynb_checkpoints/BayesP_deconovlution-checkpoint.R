suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(BayesPrism))

###----------------------------------------------
# Inputs
args<-commandArgs(TRUE)

cfrna_counts <- args[1]
ref_counts <- args[2]
ref_meta <- args[3]
biotype_key <- args[4]
output_file <- args[5]
output_file_pc <- args[6]
cores <- as.numeric(args[7])
sample_id <- args[8]

###----------------------------------------------
# Read in data
cat("-----> READING IN DATA \n")

bk.dat = read.delim(cfrna_counts,row.names=c(1),skip=1)
bk.dat <- bk.dat[,7,drop= FALSE]
colnames(bk.dat) <- sample_id

bk.dat <- bk.dat %>% t() %>% data.frame()


sc.dat <- as.matrix(fread(ref_counts),rownames=1)

sc.meta <- fread(ref_meta)
cell.type.labels = sc.meta$cell_type_collapsed

###----------------------------------------------
# Clean up reference matrix
cat("-----> CLEANING UP REFERENCE \n")
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                    species="hs", 
                                    gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                    exp.cells=5)
         
###----------------------------------------------
# Create prism                         
cat("-----> CREATE PRISM \n")                 
myPrism <- new.prism(
  reference=sc.dat.filtered, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
    outlier.fraction=0.1,
)
          
###----------------------------------------------
# Run deconvolution                          
bp.res <- run.prism(prism = myPrism, n.cores=cores)
theta <- get.fraction (bp=bp.res,
            which.theta="final",
            state.or.type="type")

###----------------------------------------------
# Save Outputs  
cat("-----> SAVING OUTPUT \n")     
output_theta = data.frame(theta)
output_theta$sample_id <- rownames(output_theta)
write.table(output_theta,output_file,sep="\t",quote=FALSE,row.names=FALSE)



###----------------------------------------------
# PROTEIN CODING ONLY 
###----------------------------------------------

###----------------------------------------------
# Extract only protein coding genes
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")

###----------------------------------------------
# Create prism                         
cat("-----> CREATE PRISM \n")                 
myPrism.pc <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
    outlier.fraction=0.1,
)
          
###----------------------------------------------
# Run deconvolution                          
bp.pc.res <- run.prism(prism = myPrism.pc, n.cores=cores)
theta.pc <- get.fraction (bp=bp.pc.res,
            which.theta="final",
            state.or.type="type")

###----------------------------------------------
# Save Outputs  
cat("-----> SAVING OUTPUT \n")    
output_theta.pc = data.frame(theta.pc)
output_theta.pc$sample_id <- rownames(output_theta.pc)                
write.table(data.frame(theta.pc),output_file_pc,sep="\t",quote=FALSE,row.names=FALSE)
    