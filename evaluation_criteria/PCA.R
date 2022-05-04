## Make PCA Plots from input expression and metadata files ##

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(gridExtra)

######################################
# DIRECTORIES - USER INPUTS REQUIRED #
######################################

# Workdir must have *metadata.csv files 
work_dir="/path/to/PCA_metadata"

# Input/Output dirs (counts_dir has the expression data, pca_dir is the output dir for PCA plots and data)

# Define main dir for corrected data 
proj2 = "/path/to/Batch_corrected_data/"

### Each batch variable has been corrected with each correction method. Select the batch-method combination below:

# What batch variable are you analyzing? Format as a string that matches the batch directory name
#batch <- 'libPrep_as_batch'
batch <- 'mission_as_batch'

# Which correction methods do you want to check? Format as a string that matches the directory name
#corrMethod = 'ComBat_seq'
#corrMethod = 'ComBat_standard'
#corrMethod = 'MBatch_AN'
#corrMethod = 'MBatch_EB'
corrMethod = 'MBatch_MP'

# Directories
counts_dir=file.path(proj2, batch, corrMethod, "Corrected_Counts")
pca_dir=file.path(proj2, batch, corrMethod, "PCA")

# Are your data log transformed? 
logt = FALSE 

#######################################
# FUNCTIONS - NO USER INPUTS REQUIRED #
#######################################

### Function to format metadata tables ###
# This function reorders the sample IDs to match the order of IDs in the expression data normCounts
# Then creates a list of factors for plotting
metFormat <- function(metTable){
  # Input: a metadata table from read.csv(), probably requires only 2 columns
  
  # Reorder sample IDs according to the expression data (named "normCounts")
  metTable2 <- as.data.frame(metTable[match(colnames(normCounts), rownames(metTable)),])
  rownames(metTable2) <- rownames(metTable)[match(colnames(normCounts), rownames(metTable))]
  colnames(metTable2) <- colnames(metTable)
  metTable <- metTable2
  
  if (dim(metTable) >= 2){
    metTableForm<-apply(metTable,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
  } else{
    metTableForm<-metTable[,1]
  }
  group_names <- paste0("(",metTableForm,")",sep = "") # human readable group names
  metTableForm <- make.names(metTableForm) # group naming compatible with R models
  names(metTableForm) <- group_names
  
  return(list('form'=metTableForm))
}

### Function to plot PCA, with given variables for color and shape ###
pcaPlot <- function(df, fileVarName, pcacolor, pcashape, metDF){
  # Inputs: df = dataframe of data to plot
  #         fileVarName = a string, the variable name to put in the file, so "GLDS" for a file called "GLDS_PCA_wlabels.png"
  #         pcacolor = a string, the metadata column name to use to color the pca points
  #         pcashape = a string, the metadata column name to use to shape the pca points
  #         metDF = the metadata dataframe 
  
  autoplot(df, data=metDF, colour=pcacolor, shape=pcashape, label=TRUE, label.size=3, size=4, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
  ggsave(file.path(pca_dir, paste0(fileVarName, '_PCA_wlabels.png')), width = 8.5, height = 6, dpi = 300)
  autoplot(df, data=metDF, colour=pcacolor, shape=pcashape, size=4, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
  ggsave(file.path(pca_dir, paste0(fileVarName, '_PCA_nolabels.png')), width = 8.5, height = 6, dpi = 300)
  
}

####################################################
# DATA AND CALCULATIONS - USER SHOULD CHECK INPUTS #
####################################################

# Import gene expression file
normCounts <- read.csv(file.path(counts_dir, "Corrected_Counts.csv"), header = TRUE, row.names = 1, stringsAsFactors = TRUE)

# Calculate PCs and write out to PCA_table.csv. Log transform if data is not already log transformed.
if (logt==FALSE) {
  exp_norm <- log2(normCounts+1)
  normCounts <- exp_norm
}

PCA_norm <- prcomp(t(normCounts), scale = FALSE)
write.csv(PCA_norm$x,file.path(pca_dir, "PCA_table.csv"), row.names = TRUE)

# Import metadata tables and format them simultaneously using metFormat() function
group <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*_cond_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
batch <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*batch_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
age <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*age_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
ar <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*animalreturn_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
dur <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*duration_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
sex <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*gender_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
lib <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*libPrep_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
miss <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*mission_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
pres <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*preservation_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
Facility <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*seqFacility_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
Parameters <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*seqParameters_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
strain <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*strain_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form
rRNA <- metFormat(read.csv(Sys.glob(file.path(work_dir,"*rRNA_metadata_Proj2.csv")), header=TRUE, row.names=1, stringsAsFactors=TRUE))$form

# Create SampleTable with imported metadata
sampleTable <- data.frame(Group=factor(group), Dataset=factor(batch), Age=factor(age), AR=factor(ar), Duration=factor(dur), 
                          Gender=factor(sex), LibPrep=factor(lib), Mission=factor(miss), Preservation=factor(pres), 
                          SeqFacility=factor(Facility), SeqParameters=factor(Parameters), Strain=factor(strain), rRNA=factor(rRNA))
# add sample names - depends on the expression and metadata files having the same IDS in the same order
rownames(sampleTable) <- colnames(normCounts)

#########################################
# MAKE PLOTS - USER SHOULD CHECK INPUTS #
#########################################

# Make PCA plots for all variables
pcaPlot(PCA_norm, 'GLDS', 'Dataset', 'Group', sampleTable)
pcaPlot(PCA_norm, 'Age', 'Age', 'Group', sampleTable)
pcaPlot(PCA_norm, 'AR', 'AR', 'Group', sampleTable)
pcaPlot(PCA_norm, 'Duration', 'Duration', 'Group', sampleTable)
pcaPlot(PCA_norm, 'Gender', 'Gender', 'Group', sampleTable)
pcaPlot(PCA_norm, 'LibPrep', 'LibPrep', 'Group', sampleTable)
pcaPlot(PCA_norm, 'Mission', 'Mission', 'Group', sampleTable)
pcaPlot(PCA_norm, 'Preservation', 'Preservation', 'Group', sampleTable)
pcaPlot(PCA_norm, 'SeqFacility', 'SeqFacility', 'Group', sampleTable)
pcaPlot(PCA_norm, 'SeqParameters', 'SeqParameters', 'Group', sampleTable)
pcaPlot(PCA_norm, 'Strain', 'Strain', 'Group', sampleTable)
pcaPlot(PCA_norm, 'rRNApct', 'rRNA', 'Group', sampleTable)

# Plot multiple PCA plots in one image
# First 6 = most important technical variables
p1 <- autoplot(PCA_norm, data=sampleTable, colour='Dataset', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p2 <- autoplot(PCA_norm, data=sampleTable, colour='Mission', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p3 <- autoplot(PCA_norm, data=sampleTable, colour='LibPrep', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p4 <- autoplot(PCA_norm, data=sampleTable, colour='SeqFacility', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p5 <- autoplot(PCA_norm, data=sampleTable, colour='Duration', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p6 <- autoplot(PCA_norm, data=sampleTable, colour='Preservation', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)

# Next 6 = biological variables, rRNA, AR and SeqParameters
p7 <- autoplot(PCA_norm, data=sampleTable, colour='Age', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p8 <- autoplot(PCA_norm, data=sampleTable, colour='Gender', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p9 <- autoplot(PCA_norm, data=sampleTable, colour='Strain', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p10 <- autoplot(PCA_norm, data=sampleTable, colour='rRNA', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p11 <- autoplot(PCA_norm, data=sampleTable, colour='AR', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)
p12 <- autoplot(PCA_norm, data=sampleTable, colour='SeqParameters', shape='Group', 
               size=2, alpha=5) + scale_shape_manual(values=c(17,16)) + theme_classic(base_size = 16)

combo1 <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)
ggsave(file.path(pca_dir, 'Combo1_PCA_nolabels.png'), plot=combo1,
       width=20, height=10, dpi=300)

combo2 <- grid.arrange(p7, p8, p9, p10, p11, p12, nrow=2)
ggsave(file.path(pca_dir, 'Combo2_PCA_nolabels.png'), plot=combo2,
       width=20, height=10, dpi=300)
