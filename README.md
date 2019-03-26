# SYNC-mutations_CRISPR
Analysis of CRISPR data from cell lines with SYNC mutations from publicly available datasets.

Data Source: DepMap, Broad Institute
    https://depmap.org/portal/download/
    
Datasets required for analysis: 
1. Genetic Dependency CRISPR (Avana) Public 19Q1 
2. Cellular Models Cell Line Metadata
3. Cellular Models Mutation Public 19Q1

CRISPR data containing batch corrected CERES score for 17,634 genes in 558 cell lines availble. CERES score is an estimate of gene's essentiality for the cell's viability. 
Score of 0 : gene is non-essential, < -0.5 : gene is essential, > 0.5 indicates KO of gene gives proliferative advantage.

Methodology & Comments

17 genes and their rare mutations identified as part of the SYNC project.  
Ref. Slides presented on 06.02.19 by Claudia Scholl. 

Following workflow applied for each gene
1. Cross referenced with mutation data of cell lines and identified cell lines containing gene mutations. (Plot 1 and 2 )
2. Distribution of CERES of gene across all cell lines evaluated. This should generally be centered around 0 and normally distributed. CERES of gene in cell lines containging non-synonymous mutations indicated. If this is below -0.5, it would imply the mutation is probably oncogenic. (Plot 3)
3. CERES data for each cell line containing any of the SYNC mutations for that gene was analysed. 
4. CERES of each gene in cell line of interest compared to the distribution of CERES of that gene across all cell lines, and pvalue computed. (Refer to Distribution_parameters.R for details of statistical test)
5. Calculation of Delta CERES = CERES of gene in cell line of interest - mean CERES of gene across cell lines (mu calculated according to Distribution_parameters.R)
6. Distribution of dCERES assessed for all genes to see if it is centered around 0. Gives an assessment of how essential the gene is. (Plot 4)
7. dCERES plotted against -log10FDR. If FDR < 0.05 and |dCERS| > 0.5,  gene is flagged as significant. Note that genes with a dCERES around -0.5 is indicative of a synthetic lethality. (Plot 5)
8. Table of genes with significantly different CERES from population mean for each cell line and common genes across all cell lines generated.
		







