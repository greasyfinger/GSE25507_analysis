# GSE25507_analysis
# Gene Expression Analysis Project

This README provides an overview of the analysis performed on a microarray dataset obtained from the Gene Expression Omnibus (GEO) with accession number GSE25507.

## Dataset Information

- **Dataset Accession Number:** [GSE25507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25507)
- **Number of Samples:** 146
  - 82 samples from individuals with autism (test group)
  - 64 samples from the control group

The dataset was loaded into an R file using the `getGEO()` function from the GEOquery library and stored in a variable named `gset`. Since only one platform (GPL570) is present, the platform index `idx` is set to 1.

## Preprocessing

In the preprocessing stage, the following steps were performed:

1. Separation of control and autism groups.
2. Setting column names for `gset` for later use in creating a top table.
3. Extraction of the microarray dataset and normalization.
4. Creation of a summary output and storage in `summ_ex.txt`.
5. Creation of a boxplot and a plot of expression values stored as `box_plt.png` and `plt_dnsty.png`, respectively.
6. Printing of `pdata` and `fdata` in an organized command line output.

## Log Transformation

Data was log-transformed to make the distribution more normal and suitable for analysis. After log transformation and re-normalization, boxplots and density plots were recalculated. The log transformation helped to make the distribution of values more centered and stabilized the variance.

## Differential Gene Expression Analysis

### Simple t-Test (Question 4)

A simple t-test was performed between the control and autism groups without using the limma package. The following steps were taken:

1. Extraction of `pdata` for control and autism groups.
2. Calculation of t-test statistics, log fold change, p-values, and adjusted p-values.
3. Generation of a volcano plot (`simple_ttest_volcano.png`) to visualize log2 fold change vs. -log10 p-values.

### Limma Package (Question 5)

Differential gene expression analysis using the limma package was performed as follows:

1. Fitting a linear model using `lmFit`.
2. Setting up contrasts of interest.
3. Calculation of statistics using `eBayes`.
4. Creation of a top table of significant genes.
5. Generation of a volcano plot (`limma_volcano.png`) to visualize the relationship between log p-value and log fold change.

## Adjusting Significance Cutoff (Question 6)

The significance cutoff was adjusted to 0.15 to include more samples for enrichment analysis. A cutoff of 0.05 is commonly chosen to maintain a 5% chance of false positives, while a cutoff of 0.15 allows a 15% chance.

## Gene Enrichment Analysis (Question 7)

Gene enrichment analysis was performed using the enrichGO function with the following steps:

1. Extraction of gene symbols from the top table.
2. Gene set enrichment analysis using the Gene Ontology Biological Processes (BP) and Homo sapiens database (org.Hs.eg.db).
3. Results were filtered based on p-value and q-value cutoffs.
4. Displaying and saving the enrichment results in `enrich_result.txt`.

Dot plots and bar plots were generated to visualize the enriched pathways, which are attached as `dot_plot.png` and `bar_plot.png`.

## Additional Information

### Top 50 Pathways (Question 8)

A subset of the top 50 pathways obtained from the enrichment analysis is provided in `top_50.txt`. These pathways mainly relate to neuronal development and growth, but also include other processes like collagen catabolic process and secondary palate development.

## Conclusion

This analysis provides insights into the gene expression differences between control and autism groups and identifies significant biological processes related to neuronal development and growth. The README summarizes the key steps and results obtained during the analysis.