![Alt text](TraianProt.png) 

<br>
<br>
TraianProt is a web-based, user-friendly proteomics data analysis platform, that enables the analysis of label free and TMT labelled data from data-dependent or data-independent acquisition mass spectrometry mode  supporting MaxQuant, MSFragger, DIA-NN, ProteoScape and Proteome Discoverer output formats. Among its functionalities  preprocessing of data, differential expression analysis, functional and protein interaction analysis are included along with the visualization of each of the aforementioned stages.


TraianProt can be accessed through the following url: https://samueldelacamara.shinyapps.io/TraianProt/

Additionally, in this repository, there is a pipeline whose files can be found in the "R" folder containing the main.R file, which drives the developed functionalities, and the files dataprocessing.R, quality_metrics.R, statisticalAnalysis.R, overviewfigures.R, enrichment.R, and Network.R, which contain the set of implemented functions. The set of files in the "Pipeline" folder serves the functions of the application in RStudio.

<br>

# Requisites
* R version >= 4.2
* Rstudio
  
- `ProteinGroups.txt` file from MaxQuant. a test dataset can be found on "data" folder.
- `combined_protein.tsv` file from MSFragger.
- `report.pg.matrix.tsv` file along with `report.tsv` file from DIA-NN in the same directory.
- Output file from Proteome Discoverer. 

<br>

# Installation

* Download the `app.R` and  `functions.R` files and load them into you R session.
* Click Run app.
* Access through the url provided above.

# Test dataset
Test datasets are available in `data` folder

<br>

# Reference:

[1] de la Cámara-Fuentes, S., Gutiérrez-Blázquez, D., Hernáez, ML., Gil, C. (2024). TraianProt: a user-friendly R shiny application for wide format proteomics data downstream analysis.  arXiv:2412.15806.
