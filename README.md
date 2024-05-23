# Triple Negative Breast Cancer Subtyping

Breast cancer is the most prevalent cancer worldwide. In 2020, 2.3 million women were diagnosed with breast cancer. Triple-negative breast cancer (TNBC) is a specific subtype of breast cancer which does not express estrogen receptor (ER), progesterone receptor (PR) and epidermal growth factor receptor 2 (HER2). TNBC shows high metastatic potential, high invasiveness, proneness to relapse and poor prognosis. Because of lack expression of ER, PR and HER2, TNBC tumors remain insensitive to classical endocrine therapy or HER2 treatment. However, TNBC stands as a heterogeneous disease that can be classified into distinct molecular subtypes using gene expression profiling leading to a subtype-specific treatment strategy. 

Here, we subtype TNBC data as 4 groups based on Jiang et al :

- LAR
- IM
- MES
- BLIS

# Prerequisites

To perform TNBC subtyping you will need :
  - R 4.1.0 with BiocManager v3.13. All packages required are listed in R_requirements.txt
  - python3 with scipy, statsmodels and pandas installed. Packages are listed in python_requirements.txt
    
# Perform subtyping

A minimal command line to run classification is :

```python3
python3 classification.py -d ${DATASET_ID} -i ${INPUT_FILE} -o ${OUTDIR} -m ${MODE}
```

Dataset ID will be used to name output files.
Mode must be "rnaseq" or "microarray".
Output folder must exists.

Other arguments are available for classification script :

- -s, --split : Number of data split to perform for crossing gene name with probes name. Depending on specs machine, do not perform split can lead to out of memory exception. Only for microarray mode. (default is 200)
- -t, --threshold : Pvalue cutoff for gene expression comparison between clusters (default is 0.05) 

# Input file

Depending on which mode is selected, expected input will differ.

## Microarray mode

For microarray mode, input file must be a four columns tab separated with no header.
First column must contain absolute path to CEL file.
Second column must contain cohorte name of the corresponding CEL file.
Third column must contain sample ID.
Fourth column must contain microarray platform experiment. Must be one of the following :

| affy_hc_g110 | affy_hg_focus | affy_hg_u133a |
| :----------: | :-----------: | :-----------: |
| affy_hg_u133a_2 | affy_hg_u133a_2 | affy_hg_u133b |
| affy_hg_u133_plus_2 | affy_hg_u95a | affy_hg_u95av2 |
| affy_hg_u95b | affy_hg_u95c | affy_hg_u95d |
| affy_hg_u95e | affy_hta_2_0 | affy_ht_hg_u133_plus_pm |
| affy_huex_1_0_st_v2 | affy_hugenefl | affy_hugene_1_0_st_v1 |
| affy_hugene_2_0_st_v1 | affy_hugene_2_1_st_v1 | affy_primeview |
| affy_u133_x3p |  |  |


For an example of microarray input file, take a look at input_microarray.txt.

## RNA-Seq mode

For RNA-Seq mode, input file must be a table of TPM count with genes in row and samples ID in columns. Gene names are expected to be rownames in the table.
Cleaning can be performed using -c option.
Genes with more than 30% of 0 among all samples count will be removed and the top 1000 genes with the highest standard deviation across samples will be conserved.

Those two parameters can be customed using -p and -n options.
Cleaning can only be performed for RNA-Seq data.
By default, cleaning is not performed.

# Output files

Output files are organised as follow :

```bash
OUTDIR
|----QC
          |NUSE.jpg (if microarray mode set)
          |RLE.jpg (if microarray mode set)

|----Clustering
          |berstein.jpg (Kmeans plot of the four groups using )

|----Boxplot
          |overExpressed.jpg
          |IL6_JAK1_STAT3.jpg
          |STATS.jpg
          |SOX.jpg
          |MES_specific.jpg
          |LAR_specific.jpg
          |FOXC1.jpg
          |AR.jpg
          |MKI67.jpg
          
|----pathwaysEnrichment
          |Cluster_1.jpg
          |Cluster_2.jpg
          |Cluster_3.jpg
          |Cluster_4.jpg

|${DATASET}_pval.csv
```

Boxplot overExpressed corresponds to specific over expressed genes for each specific subtypes.

STATS boxplot expression of STAT family genes. (STAT1, STAT2, STAT4, STAT3, STAT5A, STAT5B, STAT6)

SOX boxplot expression of SOX family genes (SOX1, SOX2, SOX3, SOX4, SOX5, SOX6, SOX7, SOX8, SOX9, SOX10, SOX11, SOX12, SOX13, SOX14, SOX15, SOX17, SOX18, SOX21, SOX30).

MES_specific boxplot expression of OGN, ADIPOQ, PLIN1, IGF1.

LAR_specific boxplot expression of ESR1, PGR, FOXA1, FOXA2, FOXA3, GATA3.

The file ${DATASET}_pval.csv contains gene expression comparison between four clusters. Each cluster is compared to the other ones.

Gathering informations leads to assign one of the four subtypes to each cluster.




