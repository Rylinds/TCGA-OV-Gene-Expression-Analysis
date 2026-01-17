# TCGA-OV Gene Expression Analysis

## Intro
The goal of this project was to compare TCGA-OV gene expression data between less severe (figo stage iia-iic) to more severe (figo stage iiia-iiib) ovarian serous cystadenocarcinoma sample groups, testing whether the more advanced-stage group exhibits genetic abnormalities indicative of disease progression.

## Methods
Cases were collected from the Genomic Data Commons (GDC) Data Portal (Heath et al., 2021). Samples were selected from the TCGA-OV project (The Cancer Genome Atlas Research Network, 2011) and limited to serous cystadenocarcinoma, nos cases with Gene Expression Quantification data produced using the STAR-Counts workflow. Only cases with open access data were included. Two groups were created based on the figo stage to define severity. Group 1 (less severe) included cases classified as stage iia-iic. Group 2 (more severe) included cases classified as stage iiia-iiib. All cases were restricted to patients with an “alive” vital status.

Sample data were downloaded from the GDC Analysis Center as CSV sample sheets, then imported in R and structured as data frames. Data retrieval and preparation were performed using the TCGAbiolinks package (Colaprico et al., 2016), which queried and downloaded the transcriptomic and clinical data. Processed transcriptomic and clinical data were stored in a structured container format designed for integrative genomic analyses using SummarizedExperiment (Morgan et al., 2023). After import, samples were screened to ensure unique subject identifiers and complete barcode information. Duplicated or overlapping samples between groups were removed, and all samples were verified to contain complete barcode information. 

After the data was imported and verified, samples were filtered using principal component analysis (PCA) to detect and remove any outliers. Genes (loci) were further filtered to remove those with excessive zero read counts across samples. Loci were excluded where more than 90% of samples in either group showed zero expression. These steps were designed to reduce data noise and improve reliability in differential expression analysis (downstream in the pipeline).

Differential gene expression between Groups 1 and 2 was calculated using the DESeq2 package (Love et al., 2014). This package uses a negative binomial model to estimate expression differences between conditions, while normalizing for sequencing depth and variability. Significance was determined using FDR adjusted p-values, with a cutoff padj < 0.05 to identify statistically significant genes. Visualization was performed in R using ggplot2 (Wickham, 2016) and ggrepel (Slowikowski, 2024) to produce PCA plots and volcano plots that summarize differential expression trends between the two sample groups.

## Results
A total of 27 samples were included before filtering, consisting of 17 Group 1 (less severe; figo stage iia-iic) and 10 Group 2 (more severe; figo stage iiia-iiib) ovarian serous cystadenocarcinoma cases. Following the first round of PCA, one Group 2 sample was identified as an outlier based on the cutoff value, PC1 > 100, and was consequently removed from further analysis (Fig. 1). The final data contained 26 samples across both groups.

Figure 1. PCA sample before removing outliers.
<img width="952" height="568" alt="Screenshot 2026-01-17 at 13 54 52" src="https://github.com/user-attachments/assets/b4207ca0-d61a-4313-a607-60354f6ac21a" />

*Each point represents an individual sample, colored by Group (Group 1 = red, less severe; Group 2 = blue, more severe). Principal components 1 and 2 explain 11.77% and 8.34% of total variance. The vertical red line marks the cutoff (PC1 = 100) used to identify an outlier sample.*

The initial dataset contained 60,660 loci, which were filtered to exclude genes with excessive zero counts across both samples, yielding 46,294 loci for differential expression analysis. Among these, 31,036 genes had non-missing false discovery rate p-adjusted values. Differential expression testing was conducted using DESeq2, comparing Group 2 (more severe) relative to Group 1 (less severe). Positive log2 (fold change) values represent genes with higher expression in Group 2, while negative values indicate downregulation in Group 2 relative to Group 1.

Using a significance threshold of padj < 0.05, a total of 17 genes were found to be differentially expressed between the two groups. 9 genes were significantly upregulated and 8 were significantly downregulated in Group 2 relative to Group 1. The top 10 most significant differentially expressed genes (ranked by p-adjusted values, smallest to largest) were: HOXA3, NTSR2, HOXA5, H19 , LINC01518, TKTL1, RIPPLY1, TRH, ORM1, ZDHHC22.

Cluster analysis with PCA revealed moderate separation between the groups along the first two principal components, which explained 10.83% and 8.91% of total variance (Fig. 2). 

Figure 2. PCA sample after removing outliers.
<img width="952" height="592" alt="Screenshot 2026-01-17 at 13 56 21" src="https://github.com/user-attachments/assets/402ce1e6-c154-4aff-b988-af61c6021c37" />

*PCA excluding the Group 2 outlier shows improved clustering. PC1 and PC2 explain 10.83% and 8.91% of variance.*

Although partial overlap remained between the groups, the analysis confirmed that sample quality and groupings were consistent after outlier removal.
A volcano plot summarizing log2 (fold change) versus -log10 (FDR p-adjusted) values showed the distribution of differentially expressed genes (Fig. 3). 

Figure 3. Distribution of log2 (fold change) and the log transformed FDR adjusted p-values across all genes.
<img width="880" height="734" alt="Screenshot 2026-01-17 at 13 57 14" src="https://github.com/user-attachments/assets/25188b83-39e1-4db3-a777-e57ab5ad77b8" />

*Each point represents a gene. The x-axis shows log₂ fold change, with genes exceeding log2 FC > 1.5 highlighted. The y-axis shows -log10(FDR-adjusted p). Triangular points indicate significantly differentially expressed genes (padj < 0.05). Blue points indicate genes exceeding both statistical and fold-change thresholds (log2 FC > 1.5 and padj < 0.05).*

Upregulated genes, such as HOXA3, NTSR2, and RIPPLY1, displayed positive fold changes. Downregulated genes, such as H19 , LINC01518, and TKTL1, displayed large negative fold changes. HOXA3 was the most significantly upregulated gene in the dataset, appearing as a distinct outlier in the volcano plot with the highest -log10 value and one of the largest positive log2. These results suggest distinct transcriptional alterations between less and more severe ovarian cancer stages, with the HOXA gene family showing strong upregulation in the more advanced disease group.

## Analysis
Among the genes most strongly upregulated in the advanced stage group (figo iiia-iiib), HOXA3 was the most transcriptionally distinct, appearing as one of the highest-magnitude positive fold changes in the volcano plot. HOXA3 belongs to the HOX family of developmental transcription factors, which regulate cell identity and differentiation during embryogenesis, but become reactivated in multiple cancers where they can promote altered epithelial phenotypes (Eoh et al., 2023). Recent research has shown that HOXA3 interacts with HOXA-AS3, a long noncoding RNA that drives epithelial-mesenchymal transition (EMT) and accelerates ovarian cancer progression (Eoh et al., 2023). Knockdown of HOXA-AS3 suppresses cell proliferation and migration while reducing EMT-related markers, suggesting that it enhances metastatic potential through regulators. The upregulation of HOXA3 in the higher stage samples of this dataset is consistent with its established role in promoting EMT-like transitions and may reflect increased invasiveness in figo stage iiic precursor lesions.

In contrast, one of the strongest downregulated genes in the advanced staged group was H19, a noncoding RNA with well-documented roles in growth regulation. Although H19 is sometimes characterized as oncogenic through mechanisms such as EMT activation and miRNA sponging (Lim et al., 2021), recent work highlights its context-dependent behavior, with multiple studies demonstrating clear tumor-suppressive functions in vivo. Mouse models lacking H19 develop larger teratocercinomas and polyp burden, showing that loss of H19 can accelerate tumor progression (Yoshimizu et al., 2008). H19 downregulation in ovarian cancer has also been associated with epithelial plasticity and disruption of imprinting patterns during malignant tissue changes (Chen et al., 2000). The decreased expression of H19 in the more severe disease group observed here aligns with prior findings that a loss of H19 correlates with more aggressive tumor phenotypes, contributing to loss of regulatory control as tumors advanced from stage ii to iii. 

Overall, the differential expression analysis identified transcriptional differences associated with disease severity and confirmed that the computational pipeline functioned as intended. The results are likely given that the major upregulated and downregulated genes correspond to known regulators of EMT, growth control, and tumor progression. However, analysis cannot establish the causal relationships of gene expression and tumor clinical severity in the data. Validation is required to establish that HOXA3 overexpression promotes serous carcinoma or whether H19 loss contributes to the shift from early to late stage ovarian cancer. Follow-up studies could include targeted knockdown or overexpression of these genes in ovarian cancer cells to clarify what regulatory mechanisms influence these transcriptional differences. Increasing the sample size could also provide more statistical strength and improve how well these findings can be generalized. 

### References
Chen, C. L., Ip, S. M., Cheng, D., Wong, L. C., & Ngan, H. Y. (2000). Loss of imprinting of the IGF-II and H19 genes in epithelial ovarian cancer. Clinical cancer research : an official journal of the American Association for Cancer Research, 6(2), 474–479.

Colaprico, A., Silva, T. C., Olsen, C., Garofano, L., Cava, C., Garolini, D., Sabedot, T., Malta, T. M., Pagnotta, S. M., Castiglioni, I., Ceccarelli, M., Bontempi, G., & Noushmehr, H. (2016). TCGAbiolinks: An R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Research, gkv1507. https://doi.org/10.1093/nar/gkv1507

Eoh, K. J., Lee, D. W., Nam, E. J., Kim, J. I., Moon, H., Kim, S. W., & Kim, Y. T. (2023). HOXA‑AS3 induces tumor progression through the epithelial‑mesenchymal transition pathway in epithelial ovarian cancer. Oncology reports, 49(3), 64. https://doi.org/10.3892/or.2023.8501

Heath, A.P., Ferretti, V., Agrawal, S. et al. The NCI Genomic Data Commons. Nature Genetics 53, 257-262 (2021). https://doi.org/10.1038/s41588-021-00791-5

Lim, Y. W. S., Xiang, X., Garg, M., Le, M. T. N., Wong, A. L.-A., Wang, L., & Goh, B.-C. (2021). The double-edged sword of H19 lncRNA: Insights into cancer therapy. Cancer Letters, 500, 253–262. https://doi.org/10.1016/j.canlet.2020.11.006

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

Morgan M, Obenchain V, Hester J, Pagès H (2023). SummarizedExperiment: SummarizedExperiment container. doi:10.18129/B9.bioc.SummarizedExperiment https://doi.org/10.18129/B9.bioc.SummarizedExperiment, R package version 1.32.0, https://bioconductor.org/packages/SummarizedExperiment 

Slowikowski, K. (2024). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.9.6, https://CRAN.R-project.org/package=ggrepel

The Cancer Genome Atlas Research Network. (2011). Integrated genomic analyses of ovarian carcinoma. Nature, 474, 609–615. https://doi.org/10.1038/nature10166
Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org

Yoshimizu, T., Miroglio, A., Ripoche, M. A., Gabory, A., Vernucci, M., Riccio, A., Colnot, S., Godard, C., Terris, B., Jammes, H., & Dandolo, L. (2008). The H19 locus acts in vivo as a tumor suppressor. Proceedings of the National Academy of Sciences of the United States of America, 105(34), 12417–12422. https://doi.org/10.1073/pnas.0801540105
