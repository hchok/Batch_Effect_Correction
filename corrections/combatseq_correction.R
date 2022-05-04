library(sva)

## Import size-normalized raw flt/gc counts, convert to matrix
proj2_flt_gc_norm_counts <- as.matrix(read.csv(Sys.glob(file.path("/path/to/Proj2_Normalized_Counts.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))

## Import metadata for raw flt/gc counts, converts to matrix
proj2_mission_metadata <- as.matrix(read.csv(Sys.glob(file.path("/path/to/all_FLT_GC_mission_metadata_Proj2.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))
proj2_libPrep_metadata <- as.matrix(read.csv(Sys.glob(file.path("/path/to/all_FLT_GC_libPrep_metadata_Proj2.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE))

## Run ComBat with different batches and covariates
combatseq_mission_flt_gc_counts <- ComBat_seq(counts=proj2_flt_gc_norm_counts, batch=proj2_mission_metadata[,1])
combatseq_libPrep_flt_gc_counts <- ComBat_seq(counts=proj2_flt_gc_norm_counts, batch=proj2_libPrep_metadata[,1])

## Set negative counts values to 0
combatseq_mission_flt_gc_counts[combat_mission_flt_gc_counts<0] <- 0
combatseq_libPrep_flt_gc_counts[combat_libPrep_flt_gc_counts<0] <- 0

## Export as .csv files
write.csv(combatseq_mission_flt_gc_counts,file.path("/path/to/mission_as_batch/ComBat_seq/Corrected_Counts/Corrected_Counts.csv"), row.names = TRUE, col.names = TRUE)
write.csv(combatseq_libPrep_flt_gc_counts,file.path("/path/to/libPrep_as_batch/ComBat_seq/Corrected_Counts/Corrected_Counts.csv"), row.names = TRUE, col.names = TRUE)



