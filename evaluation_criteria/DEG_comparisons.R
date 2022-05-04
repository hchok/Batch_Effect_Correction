# Compare FLT vs. GC DEGs in individual non-corrected GLDS datasets vs. individual corrected GLDS datasets ("within-dataset comparison")
# Compare FLT vs. GC DEGs overlaps between all pairwise GLDS datasets in noncorrected vs. corrected data ("across-dataset comparison")

## This script assumes the following file directory structure: 
## main dir
##    batch variable directory(s) (e.g. "mission_as_batch")
##        correction algorithm directory(s) (e.g. "MBatch_AN")
##        directory named "DEG_Comparisons" 

## Outputs: Writes to "DEG_Comparisons" dir in the batch variable dir: e.g. ~/mission_as_batch/DEG_Comparisons
# - total_table.csv (all DEGs for all GLDS for each batch correction method plus uncorrected)
# - match_table.csv (the counts of matching DEGs between the original data and the corrected data for each GLDS and corr method)
# - extra_table.csv (difference between match and total)
# - match_and_extra_table.csv (matching DEGs between original and corrected for each GLDS and corr method, with extra DEGs in parentheses)
# - [batchMethod]_raw_overlap_table.csv (the total after-correction overlap for all pairwise GLDS comparisons, and in parentheses the number of genes from the after-correction overlap that are preserved from the original overlap)
# - [batchMethod]_percent_overlap_table.csv (the percent of the before-correction overlap that is preserved in the after-correction overlap)
# - [batchMethod]_overlap_gridMap.pdf (the percent of the before-correction overlap that is preserved in the after-correction overlap)

## This script assumes that the original DGE tables have columns "Adj.p.value_.FLT.v.GC." and "Log2fc_.FLT.v.GC."
## and that the corrected DGE tables have columns "Adj.p.value_.FLT_47.v.GC_47." and "Log2fc_.FLT_47.v.GC_47." for all supplied GLDS numbers


library(tximport)
library(DESeq2)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(tidyHeatmap)
library(circlize)
library(S4Vectors)

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


# What batch variable are you analyzing? Format as a string that matches the batch directory name
#batch="mission_as_batch" 
batch="libPrep_as_batch"

# Which correction methods do you want to check? Format as a list of strings that match the directories in the batch variable directory
corrMethods = c('ComBat_seq', 'ComBat_standard', 'MBatch_AN', 'MBatch_EB', 'MBatch_MP')

# Where is all your batch corrected data? (This is the directory holding the batch directories, like "mission_as_batch")
batchDir <- '/path/to/Batch_corrected_data/'

# Within which GLDS dataset numbers do you want to compare FLT/GC? Format as a list of integers and strings
gldsnumbers = c(47,'48_I','48_C',137,168,173,242,'245_LAR','245_ISST')

# What LFC cutoff do you want to use? 
lfc=1 # If lfc=0, the script will just filter on pvalue<0.05

# Would you like to append a tag to the end of each output file name? Leave empty string if not. (eg "total_table.csv" could be "total_table_noLFC.csv")
fileTag = '_updatedFactors'

# Now, just run the rest of the script and look for the outputs!

#############
# READ DATA #
#############

# Read in the differential expression table for each correction method
# Assign a variable named the same as the corresponding string from the corrMethods list
# NOTE that assign() will not work within a function to assign a variable in the outer environment
for (m in corrMethods) {
  assign(m, read.csv(file.path(batchDir, batch, m, 'DGE_data_tables/differential_expression_stat.csv'), # read in the *_stat.csv table which has the updated factors and the stat columns
                     header=TRUE, row.names=1, stringsAsFactors=TRUE))
}

# Assign path to "DEG_Comparisons" output directory
outPath <- file.path(batchDir, batch, 'DEG_Comparisons')


#############
# FUNCTIONS #
#############

#######
# filterDGE()
# Function to run  pvalue and LFC filtering on input files
# Inputs:
# 1) Variable to DGE file read in as dataframe with read.csv()
# 2) boolean, whether or not it is a batch-corrected file
# 3) GLDS number (as a string), used only for corrected data
# 4,5) optional: if the adj p value and lfc columns in the uncorrected data are named differently, provide them

# Outputs:
# 2 lists: 1) adj pvalue filtered gene list, 2) adj pvalue & lfc filtered gene list
#
#######

filterDGE <- function(dgeFile, corrected=TRUE, glds) {
  
  ### Read in DGE file
  dge <- dgeFile
  
  if(corrected == FALSE) { # uncorrected data 
      # get DEGs with adj p <0.05
      p <- dge[dge$Adj.p.value_.FLT.v.GC.<0.05,]
      p <- p[is.finite(p$Adj.p.value_.FLT.v.GC.),]
      # get DEGs with adj p < 0.05, |log2FC| > 1
      pFC <- rbind(p[p$Log2fc_.FLT.v.GC.>lfc,], # up
                   p[p$Log2fc_.FLT.v.GC.< -lfc,]) # down
  }
  
  if(corrected == TRUE) { # corrected data 
    # get DEGs with adj p <0.05 
    # corrected data has specific GLDS numbers in the column names, so fill in with paste()
    p <- dge[dge[paste('Adj.p.value_.FLT_', glds, '.v.GC_', glds, '.', sep='')]<0.05,]
    p <- p[is.finite(p[[paste('Adj.p.value_.FLT_', glds, '.v.GC_', glds, '.', sep='')]]),]
    # get DEGs with adj p < 0.05, |log2FC| > 1
    pFC <- rbind(p[p[paste('Log2fc_.FLT_', glds, '.v.GC_', glds, '.', sep='')]>lfc,], # up 
                 p[p[paste('Log2fc_.FLT_', glds, '.v.GC_', glds, '.', sep='')]< -lfc,]) # down
  }
  
  return(list("Pvalue"=rownames(p), "PvalueFC"=rownames(pFC)))
  
}

#######
# compareDGE()
# Function to compare 2 input lists of genes (original and corrected)
# Inputs:
# 1) original DGE list (output from filterDGE())
# 2) corrected DGE list (output from filterDGE())
#
# Outputs:
# 3 integers: 1) original DGE count, 2) matching DGE count, 3) nonmatching DGE count, 4) list of matching genes
#######

compareDGE <- function(originalGenes, correctedGenes) {
  
  originalCount <- length(originalGenes)
  correctedCount <- length(correctedGenes)
  matching <- sum(countMatches(originalGenes, correctedGenes))
  notmatching <- correctedCount - matching
  matchingGenes <- intersect(originalGenes, correctedGenes)
  
  return(list("original"=originalCount, "match"=matching, "nonmatch"=notmatching, "matchlist"=matchingGenes))
}


#######
# RUN #
#######

colNames <- list() # make GLDS-# column names for writing files out
for (i in gldsnumbers) { colNames <- append(colNames, paste('GLDS-',i,sep=''))  }

### ONE TIME: create a table with the original (uncorrected) overlap of DEGs between GLDS datasets
orig_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(gldsnumbers) ))
colnames(orig_table) <- gldsnumbers
rownames(orig_table) <- gldsnumbers
print(orig_table)

for (g in gldsnumbers){
  for (h in gldsnumbers){
    #print(g)
    #print(h)
    # original DEG lists
    g_dge_uncorr <- filterDGE(eval(as.name((paste0('GLDS',g,sep='')))), corrected=FALSE, g)$PvalueFC
    h_dge_uncorr <- filterDGE(eval(as.name((paste0('GLDS',h,sep='')))), corrected=FALSE, h)$PvalueFC

    #overlap <- compareDGE(g_dge_uncorr, h_dge_uncorr)$match
    #print(overlap)

    #calculate and fill in before-correction overlap number
    orig_table[toString(g),][toString(h)] <- compareDGE(g_dge_uncorr, h_dge_uncorr)$match
    orig_table[toString(h),][toString(g)] <- compareDGE(g_dge_uncorr, h_dge_uncorr)$match
  }
}
print(orig_table)

write.csv(orig_table, file.path('/Users/lmsande2/Documents/Projects/BatchCorrection/Batch_Correction_Manuscript/Uncorrected_Data/Project-2_standard_correction/DGE',
                                'original_DEG_overlap_table.csv'), row.names = TRUE)
####


####
# Create and fill in total_table with total DGE counts
# Row names = the supplied correction methods with an additional row called "Original"
# Column names = the supplied GLDS dataset numbers
####
total_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(corrMethods)+1 ))
colnames(total_table) <- gldsnumbers
rownames(total_table) <- c('Uncorrected', corrMethods) 

# Fill in total counts from uncorrected datasets
for (g in gldsnumbers){
  total_table['Uncorrected',][toString(g)] <- length(filterDGE(eval(as.name((paste0('GLDS',g,sep='')))), corrected=FALSE, g)$PvalueFC)
}

# Iterate through correction methods
# For each correction method, iterate through GLDS numbers 
# Fill in total_table with total DEG counts for each correction method and each GLDS dataset
for (m in corrMethods){
  for (g in gldsnumbers){
    report <- filterDGE(eval(as.name(m)), corrected=TRUE, g)
    total_table[toString(g)][m,] <- length(report$PvalueFC)
  }
}

## Write out total_table with counts of total DEGs
colnames(total_table) <- colNames
write.csv(total_table, paste0(outPath, '/total_table', fileTag, '.csv'), row.names = TRUE)


####
# Create and fill in match_table with the counts of matching DEGs between the original data and the corrected data for each GLDS and each correction algorithm
# Create and fill in extra_table with the difference between the match number and the corrected number
# Create and fill in match_and_extra_table with original DEGs in first row, then matching values with extras in parentheses
# Row names = the supplied correction methods 
# Column names = the supplied GLDS dataset numbers
####
match_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(corrMethods)))
colnames(match_table) <- gldsnumbers
rownames(match_table) <- corrMethods

extra_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(corrMethods)))
colnames(extra_table) <- gldsnumbers
rownames(extra_table) <- corrMethods

match_and_extra_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(corrMethods)+1 ))
colnames(match_and_extra_table) <- gldsnumbers
rownames(match_and_extra_table) <- c('Uncorrected', corrMethods) 

# Fill in total DEG counts from uncorrected datasets
for (g in gldsnumbers){
  match_and_extra_table['Uncorrected',][toString(g)] <- length(filterDGE(eval(as.name((paste0('GLDS',g,sep='')))), corrected=FALSE, g)$PvalueFC)
}

# Iterate through correction methods
# For each correction method, iterate through GLDS numbers 
# Fill in match_table with the counts of matching DEGs between the original data and the corrected data
# Fill in extra_table with the difference between the match number and the corrected number
for (m in corrMethods){
  for (g in gldsnumbers){
     report <- compareDGE(
        filterDGE(eval(as.name((paste0('GLDS',g,sep='')))), corrected=FALSE, g)$PvalueFC, #original GLDS
        filterDGE(eval(as.name(m)), corrected=TRUE, g)$PvalueFC #corrected data
    )
  match_table[toString(g)][m,] <- report$match
  extra_table[toString(g)][m,] <- report$nonmatch
  match_and_extra_table[toString(g)][m,] <- paste0(report$match, ' (', report$nonmatch, ')')
  }
}

## Write out match and extra tables with counts of total DEGs
colnames(match_table) <- colNames
colnames(extra_table) <- colNames
colnames(match_and_extra_table) <- colNames
write.csv(match_table, paste0(outPath, '/match_table', fileTag, '.csv'), row.names = TRUE)
write.csv(extra_table, paste0(outPath, '/extra_table', fileTag, '.csv'), row.names = TRUE)
write.csv(match_and_extra_table, paste0(outPath, '/match_and_extra_table', fileTag, '.csv'), row.names = TRUE)

####
# Pairwise GLDS dataset comparison - how many DGEs overlap between different GLDS datasets before/after correction? "Across dataset comparison"
# For each correction method:
# Compare overlaps in DEGs before and after correction between/across GLDS datasets 
# Calculate % DEGs that overlap between each pairwise dataset comparison
####

for (m in corrMethods){
  print(m)
  # Table to hold percent overlapping DEGs preserved
  percent_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(gldsnumbers)))
  colnames(percent_table) <- gldsnumbers
  rownames(percent_table) <- gldsnumbers
  
  # Table to hold raw post-correlation DEG overlaps and preserved DEG overlaps
  raw_table <- data.frame(matrix(ncol = length(gldsnumbers), nrow = length(gldsnumbers)))
  colnames(raw_table) <- gldsnumbers
  rownames(raw_table) <- gldsnumbers
  
  for (g in gldsnumbers){
    for (h in gldsnumbers){
      # original DEG lists
      g_dge_uncorr <- filterDGE(eval(as.name((paste0('GLDS',g,sep='')))), corrected=FALSE, g)$PvalueFC
      h_dge_uncorr <- filterDGE(eval(as.name((paste0('GLDS',h,sep='')))), corrected=FALSE, h)$PvalueFC
      
      #calculate before-correction overlap number
      beforeCorr <- compareDGE(g_dge_uncorr, h_dge_uncorr)$match 
      
      # after correlation DEG lists
      g_dge_corr <- filterDGE(eval(as.name(m)), corrected=TRUE, g)$PvalueFC
      h_dge_corr <- filterDGE(eval(as.name(m)), corrected=TRUE, h)$PvalueFC
      
      # calculate after-correction overlap number
      afterCorr <- compareDGE(g_dge_corr, h_dge_corr)$match
      
      # calculate the number of DEGs that overlap after correction that are preserved from the original overlap 
      preservedNumber <- compareDGE(compareDGE(g_dge_uncorr, h_dge_uncorr)$matchlist, # uncorr overlap gene list
                                    compareDGE(g_dge_corr, h_dge_corr)$matchlist # corr overlap gene list 
                                      )$match # length of overlap between the 2 overlaps
      
      # raw_table: Save the total after-correction overlap, and in parentheses the number of genes from the after-correction overlap that are preserved from the original overlap
      # This value will be the same in g x h cell and in h x g cell 
      raw_table[toString(g),][toString(h)] <- paste0(afterCorr, ' (', preservedNumber, ')')
      raw_table[toString(h),][toString(g)] <- paste0(afterCorr, ' (', preservedNumber, ')')
      
      # percent_table: Save the percent of the before-correction overlap that is preserved in the after-correction overlap 
      # This value will be the same in g x h cell and in h x g cell
      percent_table[toString(g),][toString(h)] <- ( preservedNumber / beforeCorr )*100
      percent_table[toString(g)][toString(h),] <- ( preservedNumber / beforeCorr )*100
    }
  }
  
  percent_table[is.na(percent_table)] <- 0 # replace NaN with 0
  
  # Write out the tables to file with correction name
  write.csv(percent_table, paste0(outPath, paste0('/', toString(m), '_percent_overlap_table', fileTag, '.csv')), row.names = TRUE)
  write.csv(raw_table, paste0(outPath, paste0('/', toString(m), '_raw_overlap_table', fileTag, '.csv')), row.names = TRUE)
  
  # Make a grid plot from ComplexHeatmap and save it with correction name
  colnames(percent_table) <- colNames
  rownames(percent_table) <- colNames
  col_fun = colorRamp2(c(0,100), c("grey", "red")) # colors: the specification of c(0,100) needs to be the min and max values of the input dataframe)
  col_fun(seq(-3, 3))
  heatmap = Heatmap(as.matrix(percent_table), show_row_names=TRUE, 
                    cluster_rows = FALSE, cluster_columns = FALSE, 
                    col = col_fun,
                    rect_gp = gpar(col = "white", lwd = 2), # cell borders
                    heatmap_legend_param = list(title="Percent Overlapping\nDEGs Preserved"))
  save_pdf(heatmap, file.path(outPath,paste0(toString(m), '_overlap_gridMap', fileTag, '.pdf')),
           width = 6.8, height = 5, units = c("in", "cm", "mm"))
  
}


