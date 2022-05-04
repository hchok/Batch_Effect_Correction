# Batch Effect Correction for RNA sequencing data

This repo holds a set of scripts to run 5 batch effect correction algorithms on a combined RNA sequencing counts matrix, then run 5 evaluation criteria and score the effectiveness of the correction based on those criteria. 

**Algorithms:**
1) MBatch Empirical Bayes[1]
2) MBatch ANOVA[1]
3) MBatch Median Polish[1]
4) ComBat[2]
5) ComBat-Seq[3]

**Evaluation criteria:**
1) Batch QC
2) Principal Component Analysis (PCA)
3) Dispersion Separability Criterion (DSC)[1]
4) Differential gene expression (DGE)
5) Log Fold Change comparison (LFC)

**Gene expression data:**
(All datasets are mouse liver RNA sequencing data sourced from NASA GeneLab https://genelab.nasa.gov/)<br>
- GLDS-47
- GLDS-48
- GLDS-137
- GLDS-168
- GLDS-173
- GLDS-242
- GLDS-245

We evaluated the effectiveness of these algorithms on two technical batch variables which showed the greatest batch effect in the combined dataset: 
1) flight mission 
2) library preparation method

## Corrections:

MBatch corrections are performed in `corrections/mbatch_corrections.R`

ComBat corrections are performed in `corrections/combat_correction.R`

ComBat-Seq corrections are performed in `corrections/combatseq_correction.R`

**Input data for each script:**
1) A gene expression matrix (genes as rows and samples as columns)
2) A metadata file with 2 colums: column 1 as samples and column 2 as the batch variable for correction

**Output data:**
A batch-corrected gene expression matrix (same size as input matrix).

## Evaluation Criteria:
**BatchQC** - all BatchQC calculations are performed in the `scoring/Scoring_Evaluations.ipynb` notebook (see Scoring section for more details).

**Principal Component Analysis** - PCA plots and tables for corrected and uncorrected data are generated in the `evaluation_criteria/PCA.R` script.
- Input:<br>
			A gene expression matrix<br>
			A metadata table with variables to color the plots<br>
- Output:<br> 
      A PCA table<br> 
      PCA scatterplots colored by variables from the metadata<br>

**Dispersion Separability Criterion** - DSC plots and tables for corrected and uncorrected data are generated in the `evaluation_criteria/DSC.R` script. DSC is a method for quantifying the batch effects within a PCA-like plot by calculating the dispersion between sample centroids vs dispersion within each batch. 
- Input:<br>
			A gene expression matrix<br>
			A metadata table with variables to color the plots<br>
- Output:<br> 
      A DSC table<br> 
      DSC barplots colored by variables from the metadata<br>

**Differential Gene Expression** - DGE tables for corrected and uncorrected data are generated in the `evaluation_critera/DESeq2.R` script using the DESeq2 library. DEGs, log fold change and other metrics are calculated between the spaceflight and ground control samples for each dataset. The significant (p<0.05, |LFC|>1) DEGs are then compared between datasets in the `evaluation_critera/DEG_comparisons.R` script.
- Input:<br>
			A gene expression matrix<br>
- Output:<br> 
      A DGE table<br> 


**Log Fold Change comparison** - Pairwise dataset correlation of LFC (between spaceflight and ground control samples) plots and tables for corrected and uncorrected data are generated in the `evaluation_critera/LFC.R` script. 
- Input:<br>
			The DESeq2 tables for each dataset<br>
- Output:<br> 
      An LFC correlation table pairwise by dataset<br> 
      A heatmap plotting the values from the LFC correlation table<br>


## Scoring:

Each of the above criteria are used to quantify the effectiveness of the batch effect correction in the `scoring/Scoring_Evaluations.ipynb`. We define an effective batch effect correction algorithm as one which minimizes (or at least does not increase) the difference between samples from different batches, and which maximizes (or at least does not decrease) the expected difference between samples from different experimental biological conditions (in this case, spaceflight and ground control). 

Detailed description of each scoring mechanism is in the notebook, and inputs are the output tables from the scoring criteria scripts. Briefly: 
- skew and kurtosis are calculated similar to the BatchQC library for each sample, with the goal of maximizing the difference between pre- and post-correction. 
- PCA tables are used to calculate average distance between samples from different technical batches (minimize post-correction) or biological conditions (maximize post-correction).
- DSC tables are used to calculate the difference in DSC before and after correction for the technical variable (maximize) and the biological condition (minimize).
- LFC tables are used to calculate the average improvement in LFC correlation for pairwise datasets which come from different batches.
- DGE tables are used to calculate the average improvement in preserved DEGs for pairwise datasets which come from different batches.

The output of `scoring/Scoring_Evaluations.ipynb` is a table with scores for each evaluation criterion for each batch effect correction method and each technical variable that was evaluated.

Finally, `scoring/categorization.m` which calls `scoring/compute_exposing_hyperplanes.m`, `scoring/vert2lcon.m` and `scoring/lcon2vert.m` is used to identify the optimal correction algorithm and correction variable by geometrically probing the space of all allowable scoring functions to yield an aggregate volume-based scoring measure. 
- Input: the scoring table from `scoring/Scoring_Evaluations.ipynb`
- Output: representative scoring functions which yield optimum at one of the method-variable combinations, and a table listing volume percentages of optimum-invariant scoring function subsets. Method-variable combinations can be ranked by the volume percentage of their associated scoring function subset to identify which were most effective.


**R version:** R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"<br>
**Python version:** 3.8.8<br>
**Matlab version:** R2021a Update 3 (9.10.0.1684407) 64-bit (glnxa64)<br>


[1] MBatch: https://github.com/MD-Anderson-Bioinformatics/BatchEffectsPackage

[2] ComBat: Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics , 8(1), 118–127.

[3] ComBat-Seq: Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics, 2(3), lqaa078.

