## This script is for running DGE using the GL standard DESeq2 pipeline but staring with filtered, size-normalized, batch-corrected count data ##
## Calculates size factors on the corrected data and uses those to normalize the corrected data.

library(tximport)
library(DESeq2)
library(tidyverse)
library(STRINGdb) # for String database annotations
library(PANTHER.db) # for GOSLIM annotations


#batch="mission_as_batch"
batch="libPrep_as_batch"

#correctionMethod="MBatch_AN"
#correctionMethod="MBatch_EB"
#correctionMethod="MBatch_MP"
#correctionMethod="ComBat_standard"
correctionMethod="ComBat_seq"

# Define which organism is used in the study - this should be consistent with the name in the organisms.csv file, which matches the abbreviations used in the Panther database for each organism
organism <- "MOUSE"

## Define the location of the input data and where the output data will be printed to
# Directory containing metadata and organisms csv files
work_dir="/path/to/DGE_with_DESeq2"
# Directory containing combined, filtered and size-normalized counts table after batch-correction (batch-corrected counts data table) - actual data
counts_dir=file.path("path/to/Batch_corrected_data",
                       batch, correctionMethod, "Corrected_Counts")   


# Directory to output DGE analysis
DGE_output=file.path("/path/to/Batch_corrected_data",
                     batch, correctionMethod, "DGE_data_tables")

# Set your work_dir as your working directory
# Note: Be sure the following files are in your work_dir:
# *GLDScond_metadata*.csv
# organisms.csv


#######
# RUN #
#######
setwd(file.path(work_dir))

# Import metadata tables - note this is for Proj2, for Proj1 import the *GLDScond_metadata_Proj1.csv file
# The "I-C_ISS-LAR" file contains I/C and ISS/LAR comparisons as well as the original FLT/GC comparisons
study <- read.csv(Sys.glob(file.path(work_dir,"all_FLT_GC_GLDScond_I-C_ISS-LAR_metadata_Proj2.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)


# Format metadata tables

##### Group Formatting - GLDS-#_FLT and GLDS-#_GC
if (dim(study) >= 2){
  group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
  group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") # human readable group names
group <- make.names(group) # group naming compatible with R models
names(group) <- group_names

rm(group_names)


##### Define all the groups that will be compared, i.e. format the contrasts
contrasts <- combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))),2)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 

# Import normalized, batch-corrected counts table
## This is for the filtered, normalized, and corrected Proj2 data (be sure your corrected tables have the sample columns in the same order as the uncorrected table)
normCountsBC <- (as.matrix(read.csv(Sys.glob(file.path(counts_dir,"Corrected_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)))

# Subset the expression data to the samples in the metadata (this syntax also ensures that the expression columns are now in the same order as the rownames of the metadata)
normCountsBC <- normCountsBC[, rownames(study)]
dim(study)
dim(normCountsBC)

## Force normCounts table to only contain integers - note that although these are normalized counts, a DESeqDataset object can only be made from a count matrix composed only of integers 
normCountsBC <- ceiling(normCountsBC)

# Create SampleTable with imported metadata
sampleTable <- data.frame(condition=factor(group))
# add sample names
rownames(sampleTable) <- colnames(normCountsBC)


# make DESeqDataSet object
dds <- DESeqDataSetFromMatrix(normCountsBC, sampleTable, ~condition)
summary(dds)

##### For datasets without ERCC spike-in, redefine dds as dds_1 so that the downstream script does not have to be changed based on whether or not the dataset has ERCC spike-in
dds_1 <- dds

dds_1 <- estimateSizeFactors(dds_1)
dds_1 <- estimateDispersions(dds_1)
dds_1 <- nbinomWaldTest(dds_1)

## If the nbinomWaldTest results in rows that do not converge in beta, remove those rows as follows:
## This will have little effect because these are typically genes with very small counts and little power
dds_1 <- dds_1[which(mcols(dds_1)$betaConv),] 
summary(dds_1)

## Set up normalized counts variable as data frame from the DESeqDataSet object to be compatible with downstream commands
normCounts = as.data.frame(counts(dds_1, normalized=TRUE))

# Perform likelihood ratio test and calculate results, this is similar to ANOVA
dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

## Import organism table, which defines annotation databases for each organism
organism_table <- read.csv(file.path(work_dir,"organisms.csv"))


##### Generate annotated DGE tables

## Create annotation database for the organism used in the study
ann.dbi <- organism_table$annotations[organism_table$name == organism] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
	BiocManager::install(ann.dbi, ask = FALSE)
	library(ann.dbi, character.only=TRUE)
}

## start output tables with normalized sample expression values
reduced_output_table_1 <- normCounts 

##### Iterate through Wald Tests - this is to run pairwise comparisons of all groups. Wald Test is similar to a T-test
for (i in 1:dim(contrasts)[2]){
	res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
	res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
	colnames(res_1)<-c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("Stat_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
	reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
	rm(res_1)
}

#### Use the annotation database (built above) to add Gene Annotation columns to the output tables
# Use bioconductor annotation database to pull the following annotations
keytype = "ENSEMBL"
annot <- data.frame(rownames(reduced_output_table_1), stringsAsFactors = FALSE)
colnames(annot)[1]<-keytype
if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
	annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(reduced_output_table_1),keytype = keytype, column = "SYMBOL", multiVals = "first")
}
if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(reduced_output_table_1),keytype = keytype, column = "GENENAME", multiVals = "first")
}
if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(reduced_output_table_1),keytype = keytype, column = "ENSEMBL", multiVals = "first")
}
if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(reduced_output_table_1),keytype = keytype, column = "REFSEQ", multiVals = "first")
}
if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(reduced_output_table_1),keytype = keytype, column = "ENTREZID", multiVals = "first")
}

# Use string database to pull string annotations
string_db <- STRINGdb$new( version="11", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
string_map <- string_map[!duplicated(string_map$SYMBOL),]

# Use panther database to pull gene ontologies
annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")
pthOrganisms(PANTHER.db) <- organism
panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
panther <- na.omit(panther)
annot$GOSLIM_IDS <- panther
rm(string_db,string_map,panther,keytype)

# Calculate mean expression for each gene across all samples and add it to the output tables
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

# Calculate standard deviation for each gene across all samples and add it to the output tables
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

# Add F statistic p-value (similar to ANOVA p-value) column
# Note the F statistic p-value was calculated in the likelihood ration test performed above, we're just adding it to the output tables here
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

# Calculate mean expression and standard deviation for each gene across samples in each group and add group means and standard deviation to the output tables
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- group
group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
group_means <- group_means[-c(1),]
colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))
group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
group_stdev <- group_stdev[-c(1),]
colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))

reduced_output_table_1 <- cbind(reduced_output_table_1,group_means)

reduced_output_table_1 <- cbind(reduced_output_table_1,group_stdev)

# Remove variables that are no longer needed
rm(group_stdev,group_means,tcounts)

# add all annotation columns to the output table
reduced_output_table_1 <- cbind(annot,reduced_output_table_1)
rownames(reduced_output_table_1) <- NULL

reduced_output_table_1$GOSLIM_IDS <- vapply(reduced_output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))

write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))
write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression_stat.csv"), row.names = FALSE)
