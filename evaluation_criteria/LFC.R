# Compare LFC from DESeq2 of all genes between given GLDS datasets, either before or after batch effect correction.
# 
# User Inputs: 
# 1) Location of DGE tables "differential_expression.csv"
# 2) Whether corrected or uncorrected
#
#
#

# import libraries
library(ComplexHeatmap)
library(tidyHeatmap)

###############
# USER INPUTS #
###############

#### Define paths to DGE tables

# First, read in the non-corrected individual DGE tables: variable names must be formatted as GLDSX where X is the dataset number
# These tables have been generated using the updated DGE script with updated factors for 48 and 245 and retaining the Stat column
uncorrDir <- '/path/to/BatchCorrection/DGE_Rerun/'

GLDS47 <- read.csv(file.path(uncorrDir, 'GLDS-47/05-DESeq2_DGE/differential_expression.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS48_I <- read.csv(file.path(uncorrDir, 'GLDS-48/05-DESeq2_DGE/differential_expression_I_for_dge_comparison.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS48_C <- read.csv(file.path(uncorrDir, 'GLDS-48/05-DESeq2_DGE/differential_expression_C_for_dge_comparison.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS137 <- read.csv(file.path(uncorrDir, 'GLDS-137/05-DESeq2_DGE/differential_expression.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS168 <- read.csv(file.path(uncorrDir, 'GLDS-168/05-DESeq2_DGE/differential_expression_for_dge_comparison.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS173 <- read.csv(file.path(uncorrDir, 'GLDS-173/05-DESeq2_DGE/differential_expression.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS242 <- read.csv(file.path(uncorrDir, 'GLDS-242/05-DESeq2_DGE/differential_expression_for_dge_comparison.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS245_LAR <- read.csv(file.path(uncorrDir, 'GLDS-245/05-DESeq2_DGE/differential_expression_LAR_for_dge_comparison.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)
GLDS245_ISST <- read.csv(file.path(uncorrDir, 'GLDS-245/05-DESeq2_DGE/differential_expression_ISS_for_dge_comparison.csv'), header=TRUE, row.names=1, stringsAsFactors=TRUE)

# Which correction methods do you want to check? Format as a list of strings that match the directories in the batch variable directory
corrMethods = c('ComBat_seq', 'ComBat_standard', 'MBatch_AN', 'MBatch_EB', 'MBatch_MP')

# Where is all your batch corrected data? (This is the directory holding the batch directories, like "mission_as_batch")
batchDir <- '/path/to/Batch_corrected_data/'

# Within which GLDS dataset numbers do you want to compare FLT/GC? Format as a list of integers and strings
gldsnumbers = c(47,'48_I','48_C',137,168,173,242,'245_LAR','245_ISST')

# What batch variable are you analyzing? Format as a string that matches the batch directory name
#batch="mission_as_batch" 
batch="libPrep_as_batch"

# Where do you want the figures to output? 
# UNCORRECTED OUTPUT:
uncorrOutDir <- '/path/to/Uncorrected_Data/LFC/'

# CORRECTED OUTPUT:
# NOTE: The corrOutDir should have folders named like the corrMethods list, each with a folder inside called "LFC": eg ~/MBatch_AN/LFC
# This script will iterate through the corrMethods list and place LFC output into each of the correction algorithm dirs
corrOutDir <- file.path(batchDir, batch)


# Read in the differential expression table for each correction method
# Assign a variable named the same as the corresponding string from the corrMethods list
# NOTE that assign() will not work within a function to assign a variable in the outer environment
for (m in corrMethods) {
  assign(m, read.csv(file.path(batchDir, batch, m, 'DGE_data_tables/differential_expression_stat.csv'), 
                     header=TRUE, row.names=1, stringsAsFactors=TRUE))
}


#####
# FUNCTIONS
#####

# Get gene names from uncorrected and corrected data
uncorrGenes = row.names(eval(as.name(paste0('GLDS', gldsnumbers[1], sep='')))) # from first uncorr dataframe; will merge with others
corrGenes = row.names(eval(as.name(corrMethods[1]))) 

# Function to get a dataframe of the LFC column for a given GLDS# for corrected or uncorrected data
getLFC <- function(dge, corrected=TRUE, glds) {
  if(corrected == FALSE) lfc <- as.data.frame(eval(as.name(dge))$Log2fc_.FLT.v.GC., row.names = row.names(eval(as.name(dge)))) # uncorrected data just has 1 LFC column 
  if(corrected == TRUE) lfc <- dge[paste('Log2fc_.FLT_', glds, '.v.GC_', glds, '.', sep='')] # corrected data has pairwise LFC columns for all pairs of GLDS #s
  colnames(lfc) <- glds # rename output column to the glds integer
  return(list("LFC"=lfc))
}

# Function to create a transcriptomic LFC heatmap and a pairwise GLDS LFC correlation heatmap
lfcMatrix <- function(corrected=TRUE, batchAlg) {
  # Inputs: Required: corrected=TRUE or FALSE
  #         Optional: (if corrected=TRUE) batchAlg = the batch algorithm used, must match one of the variables of the read-in corrected data
  #         
  # Outputs: a heatmap of LFC values for all genes 
  #          a heatmap of LFC pairwise correlations between all datasets
  # If corrected=FALSE, creates a matrix of LFC from all given uncorrected GLDS datasets (must be read in as GLDS-X)
  # If corrected=TRUE, creates a matrix of the LFC from the corrected data for all given GLDS datasets, for the given batch correction algorithm

  if(corrected == FALSE){
    # create an empty matrix with rows from one of the GLDS datasets (will be merged with the rest), 0 columns 
    lfcMat <- data.frame(matrix(nrow = length(uncorrGenes), ncol = 0 ))
    rownames(lfcMat) <- uncorrGenes
    
    for (g in gldsnumbers) {
      # Fill in empty/growing matrix by merging in LFC columns from appropriate GLDS datasets
      uncorrLFC <- getLFC(paste0('GLDS',toString(g),sep=''), corrected=FALSE, g)$LFC # get LFC column
      lfcMatMerged <- merge(lfcMat, uncorrLFC,by.x = 'row.names', by.y='row.names', all.x=FALSE, all.y = FALSE) # merge into growing matrix
      rownames(lfcMatMerged) <- lfcMatMerged[,1] # reassign "row.names" column created by merge() to the row names position
      lfcMatMerged[,1] <- NULL # delete the "row.names" column 
      lfcMat <- lfcMatMerged # update growing matrix 
    }
  }
  
  if(corrected == TRUE){
    # create an empty matrix with one column per LFC vector, GLDS numbers as column names, and row names as genes 
    lfcMat <- data.frame(matrix(ncol = length(gldsnumbers), 
                                nrow = length(row.names(eval(as.name(batchAlg)))) )) # get number of genes from the corrected DGE table
    colnames(lfcMat) <- gldsnumbers
    rownames(lfcMat) <- row.names(eval(as.name(batchAlg))) # get gene names from corrected DGE table
    
    for (g in gldsnumbers) {
      lfcMat[toString(g)] <- getLFC(eval(as.name(batchAlg)), corrected=TRUE, g)$LFC
    }
  }
  
  return(list("LFCMatrix"=lfcMat))
  
}

# Function to save a complex heatmap from a matrix
heatmap <- function(m=lfcMatrix, # matrix
                    legTitle="Legend Title",
                    outFile='/path/to/fig.pdf',
                    w=6.8, # w=plot width, h=plot height
                    h=5,
                    showRowNames=FALSE,
                    clusterRowNames=FALSE,
                    clusterColNames=FALSE) {

  heatmap = Heatmap(as.matrix(m), 
                    show_row_names=showRowNames, 
                    cluster_rows=clusterRowNames, 
                    cluster_columns=clusterColNames, 
                    #col = col_fun,
                    heatmap_legend_param = list(title=legTitle))
  save_pdf(heatmap, outFile,
           width = w, height = h, units = c("in", "cm", "mm"))
}


#####
# RUN
#####

# Get a list of GLDS-# for labels
colNames <- list()
for (i in gldsnumbers) { colNames <- append(colNames, paste('GLDS-',i,sep=''))  }

##################
## UNCORRECTED  ##
## ONE TIME RUN ## 
##################

uncorrMatrix <- lfcMatrix(corrected=FALSE)$LFCMatrix
colnames(uncorrMatrix) <- colNames

#### Save LFC transcriptomic heatmap for uncorrected data -- have not updated this plot since early 2021
heatmap(
  m=uncorrMatrix,
  legTitle='LFC',
  outFile=paste0(uncorrOutDir, 'LFC_Transcriptomic_Heatmap.pdf', sep=''),
  w=6,
  h=5,
  showRowNames=FALSE,
  clusterRowNames=TRUE)

#### Save pairwise correlation heatmap for uncorrected data
heatmap(
  m=cor(uncorrMatrix, method='pearson'), # input is pairwise correlation matrix
  legTitle='Correlation',
  outFile=paste0(uncorrOutDir, 'LFC_Correlation_Heatmap.pdf', sep=''),
  w=6,
  h=5,
  showRowNames=TRUE,
  clusterRowNames=TRUE,
  clusterColNames=TRUE)

# write out
write.table(round(cor(uncorrMatrix, method='pearson'), digits=2),
            file.path(uncorrOutDir, 'cor_lfc.csv'), sep=',', quote=FALSE, row.names=TRUE, col.names=TRUE)

###############
## CORRECTED ##
###############

# Run lfcMatrix() for all correction algorithms in a for loop, save plots
for (c in corrMethods) { 
  print(c)
  lfc <- lfcMatrix(corrected=TRUE, batchAlg=c)$LFCMatrix # create LFC matrix 
  colnames(lfc) <- colNames # update column names
  print(head(lfc))
  
  corrMatrix <- cor(lfc, method='pearson')
  print(head(corrMatrix))
  
  # Save LFC transcriptomic heatmap for corrected data 
  # heatmap(
  #   m=lfc,
  #   legTitle='LFC',
  #   outFile=file.path(corrOutDir, c, '/LFC/LFC_Transcriptomic_Heatmap.pdf'),
  #   w=6,
  #   h=5,
  #   showRowNames=FALSE,
  #   clusterRowNames=FALSE)
  
  # Save pairwise correlation heatmap for corrected data
  heatmap(
    m=cor(lfc, method='pearson'), # input is pairwise correlation matrix
    legTitle='Correlation',
    outFile=file.path(corrOutDir, c, '/LFC/LFC_Correlation_Heatmap.pdf'),
    w=6,
    h=5,
    showRowNames=TRUE,
    clusterRowNames=TRUE,
    clusterColNames=TRUE)
  
  # write out
  write.table(round(corrMatrix, digits=2), file.path(corrOutDir, c, '/LFC/cor_lfc.csv'), sep=',', quote=FALSE, row.names=TRUE, col.names=TRUE)

  
  }




