library(MBatch)


# Location of input data
normdata <- "path/to/Combined_SizeNormalized_Counts"

# Location of metadata files 
meta_dir="/path/to/metadata"

# Batch (used to assign output dir)
#batchCorrected = "libPrep"
batchCorrected = "mission"

# Read in expression file as matrix
expMat <- as.matrix(read.delim("/path/to/expression/csv", sep=",", header=TRUE, row.names=1, stringsAsFactors=TRUE))

# Read in metadata file with sample names as the first column, rename first column to "Sample"
meta <- read.delim("/path/to/meta/csv", sep=",",header=TRUE, stringsAsFactors=TRUE)
names(meta)[1] <- 'Sample'


# Create a mbatchData object
mbatchData <- mbatchLoadStructures(theGeneMatrix = expMat, theBatchDataframe = meta)


# Apply mBatch corrections

## MBatch - EB
eb <- EB_withParametricPriors(theBeaData=mbatchData,
                                         theBatchIdsNotToCorrect=c(""), 
                                         theDoCheckPlotsFlag=FALSE,
                                         theBatchType="batch",
                                         theThreads=1,
                                         theWriteToFile=FALSE)

# Convert negatives to 0
eb[eb<0] <- 0


## MBatch - AN
an <- AN_Adjusted(theBeaData=mbatchData,
                             theBatchType="batch",
                             theWriteToFile=FALSE)

# Convert negatives and nan to 0 
an[an<0] <- 0
an[is.nan(an)] <- 0


## MBatch - MP
mp <- MP_ByBatch(theBeaData=mbatchData,
                       theBatchType="batch",
                       theWriteToFile=FALSE)

# Convert negatives and nan to 0 and add 1 pseudocount
mp[mp<0] <- 0
mp[is.nan(mp)] <- 0
mp <- mp + 1

write_dir <- "/path/to/Batch_corrected_data/"


write.csv(eb,file.path(write_dir, paste(batchCorrected, "as_batch", sep="_"), "MBatch_EB/Corrected_Counts", "Corrected_Counts.csv"), row.names = TRUE)
write.csv(an,file.path(write_dir, paste(batchCorrected, "as_batch", sep="_"), "MBatch_AN/Corrected_Counts", "Corrected_Counts.csv"), row.names = TRUE)
write.csv(mp,file.path(write_dir, paste(batchCorrected, "as_batch", sep="_"), "MBatch_MP/Corrected_Counts", "Corrected_Counts.csv"), row.names = TRUE)



