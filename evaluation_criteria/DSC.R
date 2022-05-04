# Calculate MBatch dispersion separability criterion (DSC) for a given gene expression dataset and set of batches.
# 
# Inputs: 
# - gene expression file
# - batches file 
# - output path

library(MBatch)

####
# USER INPUTS
####

###### Path to expression matrix (CSV)

# Define main dir for corrected data to shorten file paths
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

pathToExp <- file.path(proj2, batch, corrMethod, "Corrected_Counts/Corrected_Counts.csv")
theOutputDir <- file.path(proj2, batch, corrMethod, "DSC")


###### Path to metadata file (CSV) 
# Requires at least 2 columns, "Sample" and a batch, can have multiple batch columns
pathToBatches <- "/path/to/all_metadata_Proj2_updatedFactors.csv"


####
# END USER INPUTS - RUN SCRIPT
####

# Read in data
# Expression
expression <- read.delim(pathToExp, sep=",",
                         header=TRUE, row.names=1, stringsAsFactors=TRUE)
dim(expression)

# Metadata 
batches <- read.delim(pathToBatches, sep=",",
                          header=TRUE, stringsAsFactors=TRUE)
head(batches)

# Create BEA_DATA object (using empty df for covariates)
myData <- new("BEA_DATA", as.matrix(expression), batches, data.frame())
#myData


# PCA (gives DSC value)
# DSC values for top 4 PCs are in the output dir in files like:
# /theOutputDir/Batch/ManyToMany/ANY_Comp1_Comp2_DSC.txt
PCA_Regular_Structures(theData=myData,
                       theTitle='PCA',
                       theOutputPath=theOutputDir,
                       theBatchTypeAndValuePairsToRemove=NULL,
                       theBatchTypeAndValuePairsToKeep=NULL,
                       theListOfComponentsToPlot=c(1, 2),
                       theDoDSCFlag=TRUE,
                       theDoDscPermsFileFlag=TRUE,
                       theDSCPermutations=1000,
                       theDSCThreads=1,
                       theMinBatchSize=2,
                       theJavaParameters="-Xms2000m",
                       theSeed=runif(1), # random number generator
                       theMaxGeneCount=20000)

# DF to hold DSC values for each batch
dscDF <- data.frame(matrix(ncol = length(colnames(batches[-1])), nrow = 1))
colnames(dscDF) <- colnames(batches)[-1]
rownames(dscDF) <- "Overall DSC"

for (b in colnames(batches)[-1]){

  # Read in PCA annotations file to get DSC for each batch
  annot <- read.delim(file.path(theOutputDir, b, "ManyToMany/PCAAnnotations.tsv"))
  
  dscDF[toString(b)]["Overall DSC",] <- annot[3,][4]
}

write.csv(dscDF, file.path(theOutputDir, 'DSC_table.csv'), row.names = TRUE)

## Source for factors idea: https://stackoverflow.com/questions/32445106/colouring-barplot-in-r-according-to-variable-of-data-set

# make plotting df 
plotDF <- as.data.frame(t(dscDF))
plotDF$`Overall DSC` <- as.numeric(plotDF$`Overall DSC`)
plotDF <- format(plotDF, digits=1) # round
plotDF$`Overall DSC` <- as.numeric(plotDF$`Overall DSC`)

colors <- c("darkseagreen", "darksalmon")
factors <- factor(c("biological", #age
                    "technical", #animalreturn
                    "technical", #dataset
                    "biological", #condition
                    "biological", #duration
                    "biological", #gender
                    "technical", #libPrep
                    "technical", #mission
                    "technical", #preservation
                    "technical", #seqFacility
                    "technical", #seqParameters
                    "biological")) #strain
pdf(file.path(theOutputDir,"DSC_barplot.pdf")) ; par(mar=c(7,4,4,1)+.1) # margins order: 'bottom', 'left', 'top', 'right'.
b <- barplot(plotDF$'Overall DSC', names.arg=row.names(plotDF), main='',
        col=as.vector(colors[factors]), ylim=c(0,1),
        #ylim=c(0,max(plotDF$'Overall DSC')), 
        ylab="Overall DSC", las=2)
text(x=b, y=plotDF$`Overall DSC`, label=plotDF$`Overall DSC`, pos=3, 
     cex=0.8, col="black")

dev.off()


