![Alt text](TraianProt.png) 

In this repository, the downstreaming analysis code for the quantitative proteomics analysis tool can be found.
<br>
<br>
TraianProt is a web-based, user-friendly proteomics data analysis platform, that enables the analysis of label free and TMT labelled data from data-dependent or data-independent acquisition mass spectrometry mode  supporting MaxQuant, MSFragger, DIA-NN and Proteome Discoverer output formats. Its functionalities include preprocessing of data, differential expression analysis, functional and protein interaction analysis, along with visualization of each of the aforementioned stages.

Additionally, in this repository, there is a pipeline whose files can be found in the "R" folder containing the main.R file, which drives the developed functionalities, and the files dataprocessing.R, quality_metrics.R, statisticalAnalysis.R, overviewfigures.R, enrichment.R, and Network.R, which contain the set of implemented functions. The set of files in the "Pipeline" folder serves the functions of the application in RStudio.

<br>

# Requisites
* R version >= 4.2
* Rstudio
  
- "ProteinGroups.txt" file from MaxQuant. a test dataset can be found on "data" folder.
- "combined_protein.tsv" file from MSFragger.
- "report.pg.matrix.tsv" file along with "report.tsv" file from DIA-NN in the same directory.
- Output file from Proteome Discoverer. 

<br>

# Outputs
1. Quality metrics:
- Boxplots
- Dispersion plot
- Scatter plots
- Pre/post imputation plots
- Histograms
- Principal Component Analysis
- Correlation plot
- Q-Q plot
2. Overview figures:
- Volcano plot
- Heatmap
- Protein Intensity plot
3. Enrichment:
- Dotplot
- Barplot
- Manhattan plot
4. Interactions:
- Interaction networks
5. Data:
- *data set* after preprocessing
- *data set* after statistical analysis
- *data set* after functional analsysis.
  
# Credits
Packages that make it work:  gprofiler2, igraph, limma, clusterprofiler, enrichplot, stringdb.
