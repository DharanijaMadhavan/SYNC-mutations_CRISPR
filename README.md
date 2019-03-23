# SYNC-mutations_CRISPR
Analysis of CRISPR data from cell lines with SYNC mutations

Data Source: DepMap, Broad Institute
    https://depmap.org/portal/download/
    
Datasets required for analysis: Genetic Dependency CRISPR (Avana) Public 19Q1, Cellular Models Cell Line Metadata,Cellular Models Mutation Public 19Q1

CRISPR data containing batch corrected CERES score for 17,634 genes in 558 cell lines availble. CERES score is an estimate of gene's essentiality for the cell's viability. 
Score of 0 : Gene is non-essential, < -0.5 : Gene is essential, > 0.5 indicates KO of Gene gives proliferative advantage.

Methodology & Comments
17 genes and their rare mutations identified as part of the SYNC project. Ref. Slides presented on 06.02.19 by Claudia Scholl. Following workflow applied for each gene
1. Cross referenced with mutation data of cell lines and identified cell lines containing gene mutations. (Plot 1 and 2 )
2. Distribution of CERES of gene across all cell lines evaluated. This should be centered around 0 and normally distributed. CERES of gene in cell lines containging non-synonymous mutations highlighted. If this is below -0.5, it would imply the mutation is probably oncogenic. (Plot 3)
3. CERES data for each cell line containing any one of the SYNC mutations was analysed. 
4. CERES of each gene in cell line of interest compared the distribution of CERES of that gene across all cell lines. pvalue calculated. (Refer to Distribution_parameters.R for details of statistical test)
5. Delta CERES = CERES of gene in cell line of interest - mean CERES of gene across cell lines
6. Distribution of dCERES assessed for all genes to see if it is centered around 0. (Plot 4)
7. dCERES plotted against -log10FDR. If FDR < 0.1 gene is flagged as significant. Note that this combined with a dCERES around -0.5 is indicative of a synthetic lethality. (Plot 5)
8. Table of synthetically lethal genes for each cell line and common genes across all cell lines generated.
		







