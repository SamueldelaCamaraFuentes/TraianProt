if(!require(shiny)){install.packages("shiny")}
library(shiny)
if(!require(shinyWidgets)){install.packages("shinyWidgets")}
library(shinyWidgets)
if(!require(shinydashboard)){install.packages("shinydashboard")}
library(shinydashboard)
if(!require(shinyjs)){install.packages("shinyjs")}
library(shinyjs)
if(!require(DT)){install.packages("DT")}
library(DT)
if(!require(BiocManager)){install.packages("BiocManager")}
library(BiocManager)
options(repos = BiocManager::repositories())
setwd("E:/SAMUEL/Mi app")
source("functions.R")
options(shiny.maxRequestSize = 30*1024^2)



ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "TraianProt"),
  sidebar <- dashboardSidebar(
    
    #Home
    conditionalPanel(condition = "input.main_tabs == 'Home_introduction'",
      sidebarMenu(menuItem('Home', icon=icon("home"), selected = TRUE, tabName = "home"))),
    
    #Data handling
    conditionalPanel(condition = "input.main_tabs == 'data_handling'",
                     sidebarMenu(
                       menuItem("File Input", tabName = "file input", icon = icon("upload"),
                                fileInput(inputId = "file", 
                                          label = h3("File input"),
                                          multiple = FALSE, 
                                          accept = c("text/csv",
                                                     ".csv")),
                                selectInput(inputId = "comptplatform", 
                                            label = "Platform",
                                            choices = list("MaxQuant" = 1, 
                                                           "MSFragger" = 2,
                                                           "DIA-NN" = 3,
                                                           "Proteome Discoverer" = 4,
                                                           "ProteoScape" = 5),
                                            selected = 4),
                                selectInput(inputId = "labeltype", 
                                            label = "Type",
                                            choices = list("Label free" = 1, 
                                                           "TMT" = 2,
                                                           "SILAC" = 3),
                                            selected = 1),
                                selectInput(inputId = "organismfocus", 
                                            label = "Organism",
                                            choices = list("Candida albicans" = 1, 
                                                           "Other" = 2),
                                            selected = 2),
                                radioButtons(inputId = "intensitymode",
                                             label = "Choose quantification",
                                             choices = c("Intensity" = "int",
                                                         "LFQ intensity" = "lfq",
                                                         "Spectral Count" = "spcount"),
                                             selected = "lfq"),
                                
                                textInput(inputId = "condition1",
                                          label = "Control condition",
                                          value = "Regexp condition 1"),
                                
                                textInput(inputId = "condition2",
                                          label = "Treatment condition",
                                          value = "Regexp condition 2"),
                                
                                numericInput(inputId = "repcond1",
                                             label = "Control samples",
                                             value = 4),
                                
                                numericInput(inputId = "repcond2",
                                             label = "Treatment samples",
                                             value = 4)),
                       menuItem("Pre-processing", icon = icon("th"), tabName = "preprocessing",
                                h4("Filtering by samples"),
                                numericInput(inputId = "minfiltro",
                                             label = "Minimum values for filtering per condition",
                                             value = 0),
                                h4("Additional filtering step"),
                                numericInput(inputId = "numberuniquepep",
                                             label = "Minimum number of unique peptides",
                                             value = 1),
                                switchInput(inputId = "uniquefilter",
                                            label = "Unique peptides filter",
                                            value = FALSE,
                                            onLabel = "Yes",
                                            offLabel = "No"),
                                
                                h3("Normalization"),
                                switchInput(inputId = "normchoice",
                                            label = "Normalization",
                                            value = FALSE,
                                            onLabel = "Yes",
                                            offLabel = "No"),
                                textInput(inputId = "normoptions",
                                          label = "Normalization method",
                                          value = "Median centering"),
                                helpText("Note: mean, median, trimMean and vsc"),
                                
                                h3("Imputation"),
                                selectInput(inputId = "imputation", 
                                            label = "Imputation method",
                                            choices = list("No imputation" = 3,
                                                           "K Nearest Neighbors" = 1, 
                                                           "Normal distribution" = 2),
                                            selected = 3),
                                actionButton( inputId = "displaytable",
                                              label = "Display table"),
                                selectInput(inputId = "displayproteins", 
                                            label = "Choose a protein set to display",
                                            choices = list("Common proteins between conditions" = 1,
                                                           "Exclusive control proteins" = 2,
                                                           "Exclusive treatment proteins" = 3),
                                            selected = 1)),
                       menuItem("Venn Diagram", icon = icon("th"), tabName = "venn",
                                strong(h4("Condition 1")),
                                textInput(inputId = "condition1venn",
                                          label = "Control condition:",
                                          value = "Control"),
                                textInput(inputId = "color1",
                                          label = "Color:",
                                          value = "blue"),
                                strong(h4("Condition 2")),
                                textInput(inputId = "condition2venn",
                                          label = "Treatment condition:",
                                          value = "Treatment"),
                                textInput(inputId = "color2",
                                          label = "Color:",
                                          value = "maroon"),
                                selectInput(inputId = "vennextension",
                                            label = "File type:",
                                            choices = c("tiff", "jpeg", "png"),
                                            selected = "pdf"),
                                selectInput(inputId = "vennquality",
                                            label = "Quality:",
                                            choices = c("High" = "retina",
                                                        "Medium" = "print",
                                                        "low" = "screen"),
                                            selected = ""),
                                textInput(inputId = "venn_name",
                                          label = "Filename Venn",
                                          value = "Venn"),
                                h5("Venn Diagram download"),
                                downloadButton(outputId = "downloadvenn",
                                               label = "Download Venn",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;")),
                       
                       menuItem("Download tables", icon = icon("download"), tabName = "Download",
                                textInput(inputId = "filenamedownloadcomm",
                                          label = "Filename proteins",
                                          value = "data"),
                                downloadButton(outputId = "downloaddfcomm", icon("download"),
                                               label = "Download",
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                textInput(inputId = "filenamedownloadcontroluniques",
                                          label = "Filename unique control proteins",
                                          value = "data"),
                                downloadButton(outputId = "downloaddfcontroluniques", icon("download"),
                                               label = "Download",
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                textInput(inputId = "filenamedownloadtreatmentuniques",
                                          label = "Filename treatment unique proteins",
                                          value = "data"),
                                downloadButton(outputId = "downloaddftreatmentuniques", icon("download"),
                                               label = "Download",
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"))
                     )),
    
    
    #Quality metrics
    conditionalPanel(condition = "input.main_tabs == 'quality'",
                     sidebarMenu(
                       menuItem( "Quality metrics", icon = icon("chart-line"),
                                 strong(h4("Distribution plots")),
                                 selectInput(inputId = "displaydistplots", 
                                             label = "Select a plot",
                                             choices = list("Boxplot" = 1, 
                                                           "Dispersion plot" = 2),
                                                            selected = 1),
                                 strong(h4("Imputation plot")),
                                 selectInput(inputId = "displayimpplots", 
                                             label = "Select a plot",
                                             choices = list("Pre imputation" = 1, 
                                                            "Post imputation" = 2),
                                             selected = 1),
                                 strong(h4("Normality plots")),
                                 textInput(inputId = "histsample",
                                           label = "Sample for histogram:",
                                           value = "LOG2.WT1" ),
                                 textInput(inputId = "histcol",
                                           label = "Color for histogram:",
                                           value = "maroon"),
                                 textInput(inputId = "histtitle",
                                           label = "Title for histogram:",
                                           value = "Histogram of intensities" ),
                                 textInput(inputId = "qqcolor",
                                           label = "Color for QQplot:",
                                           value = "blue"),
                                 selectInput(inputId = "displaynormplots", 
                                             label = "Select a plot",
                                             choices = list("Histogram" = 1, 
                                                            "Q-Q plot" = 2),
                                             selected = 1),
                                 strong(h4("Dimension reduction")),
                                 selectInput(inputId = "dimreduction", 
                                             label = "Select a plot",
                                             choices = list("PCA" = 1, 
                                                            "t-SNE" = 2),
                                             selected = 1),
                                 numericInput(inputId = "perplexitytsne",
                                              label = "perplexity tsne",
                                              value = 1),
                                 strong(h4("Correlation plots")),
                                 textInput(inputId = "scatsample1",
                                           label = "Sample1 for scatter:",
                                           value = "LOG2.WT1" ),
                                 textInput(inputId = "scatsample2",
                                           label = "Sample2 for scatter:",
                                           value = "LOG2.WT2"),
                                 textInput(inputId = "dispmethod",
                                           label = "Display method:",
                                           value = "shade"),
                                 selectInput(inputId = "displaycorrplots", 
                                             label = "Select a plot",
                                             choices = list("Scatter plot" = 1, 
                                                            "Correlation plot" = 2),
                                             selected = 1)),
                       menuItem("Download plots", icon = icon("download"), tabName = "Download",
                       strong(h4("Quality metrics download")),
                       selectInput(inputId = "qualitymetricsplotchoice", 
                                   label = "Select a plot",
                                   choices = list("Boxplot" = 1, 
                                                  "Dispersion plot" = 2,
                                                  "Pre imputation" = 3,
                                                  "Post imputation" = 4,
                                                  "Histogram" = 5,
                                                  "Q-Q plot" = 6,
                                                  "Scatter plot" = 7,
                                                  "Correlation plot" = 8,
                                                  "PCA" = 9,
                                                  "t-SNE" = 10),
                                    selected = 1),
                       textInput(inputId = "qualitymetricsfilename",
                                 label = "Filename",
                                 value = "quality metrics plot"),
                       downloadButton(outputId = "downloadqmplot", icon("download"),
                                      label = "Download",
                                      style="display: block; margin: 0 auto; width: 200px; color:black;"))
                     )),
    
    #Diferential analysis 
    conditionalPanel(condition = "input.main_tabs == 'statistics'",
                     sidebarMenu(
                       
                       menuItem( "Differential analysis", icon = icon("table"),
                                 h3("Tipe of test"),
                                 selectInput(inputId = "displaytest", 
                                             label = "Choose a test to analyze:",
                                             choices = list("Simple t test approach" = 1, 
                                                            "Limma approach" = 2,
                                                            "Wilcoxon test" = 3),
                                             selected = 2),
                                 switchInput(inputId = "PSMaware",
                                             label = "PSMs correction",
                                             value = TRUE,
                                             onLabel = "Yes",
                                             offLabel = "No"),
                                 selectInput(inputId = "proteins", 
                                             label = "Choose a way to proceed:",
                                             choices = list("Common proteins" = 1, 
                                                            "All protein" = 2),
                                             selected = 2),
                                 
                                 switchInput(inputId = "testid",
                                             label = "Test",
                                             value = FALSE,
                                             onLabel = "Paired t-test",
                                             offLabel = "Two sample t-test"),
                                 h3("Choice:"),
                                 selectInput(inputId = "statselected", 
                                             label = "Select a statistic:",
                                             choices = list("P value" = 1, 
                                                            "Q value" = 2),
                                             selected = 1),
                                 numericInput(inputId = "LogFCup",
                                              label = "Log2FC threshold upregulated:",
                                              value = 0.585),
                                 numericInput(inputId = "LogFCdown",
                                              label = "Log2FC threshold downregulated:",
                                              value = -0.585),
                                 numericInput(inputId = "sigcutoff",
                                              label = "Statistic threshold:",
                                              value = 0.05),
                                 selectInput(inputId = "pvaladj", 
                                             label = "Choose a p-value adjustment:",
                                             choices = list("FDR" = "fdr",
                                                            "Bonferroni" = "bonferroni",
                                                            "BH" = "BH",
                                                            "Hochberg" = "hochberg"),
                                             selected = "fdr")),
                                 
                       menuItem("Download tables", icon = icon("download"), tabName = "Download",
                                textInput(inputId = "namedownload",
                                          label = "Filename",
                                          value = "differential data"),
                                downloadButton(outputId = "downloaddif", icon("download"),
                                               label = "Download",
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"))
                       
                     )),
    
    #Differential analysis plots 
    conditionalPanel(condition = "input.main_tabs == 'diffplots'",
                     sidebarMenu(
                       menuItem( "Differential analysis plots", icon = icon("chart-area"),
                                 strong(h4("Volcano plot")),
                                 numericInput(inputId = "labelprots",
                                              label = "Point size",
                                              value = 10),
                                 textInput(inputId = "volcanotitle",
                                           label = "Insert a title",
                                           value = "Treatment vs Control"),
                                 strong(h4("Heatmap")),
                                 textInput(inputId = "heatmaptitle",
                                           label = "Insert a title",
                                           value = "Treatment vs Control"),
                                 strong(h4("Protein Intensity")),
                                 textInput(inputId = "proteinname",
                                           label = "Insert a title",
                                           value = "Treatment vs Control")
                                 ),
                       
                       menuItem("Download plot options", icon =  icon("download"),
                                selectInput(inputId = "difextension",
                                            label = "File type:",
                                            choices = c("tiff", "pdf", "jpeg", "png"),
                                            selected = "pdf"),
                                selectInput(inputId = "difquality",
                                            label = "Quality:",
                                            choices = c("High" = "retina",
                                                        "Medium" = "print",
                                                        "low" = "screen"),
                                            selected = ""),
                                textInput(inputId = "name_download_volcano",
                                          label = "Filename Volcano",
                                          value = "Volcano"),
                                textInput(inputId = "name_download_heatmap",
                                          label = "Filename Heatmap",
                                          value = "Heatmap"),
                                h5("Volcano download"),
                                downloadButton(outputId = "downloadvolcano",
                                               label = "Download volcano",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                h5("Heatmap download"),
                                downloadButton(outputId = "downloadheatmap",
                                               label = "Download heatmap",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                h5("Differential Heatmap download"),
                                downloadButton(outputId = "downloaddifheatmap",
                                               label = "Download heatmap",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;")
                       )
                     )),
    
    #Functional analysis
    conditionalPanel(condition = "input.main_tabs == 'fanalisis'",
                     sidebarMenu(
                       menuItem( "Functional analysis", icon = icon("table"),
                                 textInput(inputId = "organism",
                                           label = "Organism:",
                                           value = "hsapiens"),
                                 textInput(inputId = "proteinid",
                                           label = "Type of Protein ID:",
                                           value = "ENSG"),
                                 numericInput(inputId = "userthreshold",
                                              label = "enrichment threshold:",
                                              value = 0.05),
                                 switchInput(inputId = "backgroundset",
                                             label = "Background",
                                             value = FALSE,
                                             onLabel = "Yes",
                                             offLabel = "No")
                       ),
                       menuItem( "Plots", icon = icon("chart-area"),
                                 strong(h4("Dotplot")),
                                 textInput(inputId = "dotplottitle",
                                           label = "Title:",
                                           value = "Treatment vs Control"),
                                 numericInput(inputId = "categorydotplot",
                                              label = "Number of terms to display:",
                                              value = 10),
                                 numericInput(inputId = "fontsizedotplot",
                                              label = "Font size:",
                                              value = 16),
                                 strong(h4("Barplot")),
                                 textInput(inputId = "barplottitle",
                                           label = "Title:",
                                           value = "Treatment vs Control"),
                                 numericInput(inputId = "categorybarplot",
                                              label = "Number of terms to display:",
                                              value = 10),
                                 numericInput(inputId = "fontsizebarplot",
                                              label = "Font size:",
                                              value = 10),
                                 strong(h4("Manhattan plot")),
                                 actionButton( inputId = "manhattan",
                                               label = "Render plot")
                       ),
                       
                       menuItem("Download plot options", icon = icon("download"),
                                selectInput(inputId = "funcextension",
                                            label = "File type:",
                                            choices = c("tiff", "pdf", "jpeg", "png"),
                                            selected = "pdf"),
                                selectInput(inputId = "funcquality",
                                            label = "Quality:",
                                            choices = c("High" = "retina",
                                                        "Medium" = "print",
                                                        "low" = "screen"),
                                            selected = ""),
                                textInput(inputId = "name_download_func",
                                          label = "Filename",
                                          value = "functional analysis  data"),
                                h5("Dotplot download"),
                                downloadButton(outputId = "downloaddotplot",
                                               label = "Download",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                h5("Barplot download"),
                                downloadButton(outputId = "downloadbarplot",
                                               label = "Download",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;")),

                       
                       menuItem("Download tables", icon = icon("download"), tabName = "Download",
                                textInput(inputId = "name_download",
                                          label = "Filename",
                                          value = "functional analysis  data"),
                                downloadButton(outputId = "downloadfunc", icon("download"),
                                               label = "Download",
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"))
                     )),
    
    #Interactions 
    conditionalPanel(condition = "input.main_tabs == 'interactions'",
                     sidebarMenu(
                       
                       menuItem("STRINGdb", icon = icon("atom"),
                                strong(h4("Species taxon id")),
                                numericInput(inputId = "taxonid",
                                             label = "Taxon id",
                                             value = 237561),
                                strong(h4("Score threshold")),
                                numericInput(inputId = "scthreshold",
                                             label = "Score threshold",
                                             value = 400),
                                strong(h4("Up-regulated network")),
                                actionButton( inputId = "upnet",
                                              label = "Render plot"),
                                strong(h4("Down-regulated network")),
                                actionButton( inputId = "downnet",
                                              label = "Render plot")),
                       
                       menuItem("Download plot options", icon =  icon("download"),
                                selectInput(inputId = "intextension",
                                            label = "File type:",
                                            choices = c("tiff", "pdf", "jpeg", "png"),
                                            selected = "pdf"),
                                selectInput(inputId = "intquality",
                                            label = "Quality:",
                                            choices = c("High" = "retina",
                                                        "Medium" = "print",
                                                        "low" = "screen"),
                                            selected = ""),
                                textInput(inputId = "name_download_int",
                                          label = "Filename",
                                          value = "functional analysis  data"),
                                h5("Up-regulated download"),
                                downloadButton(outputId = "downloadupnetwork",
                                               label = "Download",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                h5("Down-regulated download"),
                                downloadButton(outputId = "downloaddownnetwork",
                                               label = "Download",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"),
                                downloadButton(outputId = "downloadreport",
                                               label = "Download",
                                               icon = icon("download"),
                                               style="display: block; margin: 0 auto; width: 200px; color:black;"))
                     ))
  ),
  
  
  dashboardBody(
    tabsetPanel(id = "main_tabs",
                
                #Home tab
                tabPanel(
                  title = "Home",
                  div(
                    imageOutput("home_img", inline = TRUE),
                    style = "text-align: center; margin-left: 5%;"  # Adjust the margin as needed
                  ),
                  div(
                    style = "
                    border-radius: 10px;
                    background-color: #f9f9f9;
                    border: 2px solid #ddd;
                    padding: 20px;
                    margin-top: 20px;
                    width: 50%;
                    box-shadow: 0px 4px 10px rgba(0,0,0,0.1);
                    text-align: center;
                    margin-left: auto;
                    margin-right: auto;
                    font-family: 'Arial', sans-serif;
                  ", 
                    # Text Styling
                    h4("Welcome to TraianProt", style = "color: #2c3e50; font-weight: bold;"),
                    p("TraianProt is a web-based, user-friendly proteomics data analysis platform, that enables the analysis of label free and TMT labelled data 
                      from data-dependent or data-independent acquisition mass spec-trometry mode  supporting MaxQuant, MSFragger, DIA-NN, ProteoScape and Proteome Discoverer output formats.", 
                      style = "color: #34495e; font-size: 16px;"
                    )
                  ),
                  value = "Home_introduction",
                  icon = icon("home")
                ),

                
                #Data handling tab
                
                tabPanel(title = "Preprocessing",
                         value = "data_handling",
                         icon = icon("table"),
                         
                         column(width = 12, DT::dataTableOutput("file")),
                         fluidPage(
                           fluidRow(
                             column(3, h3("Columns"),tableOutput("LOG2.names")),
                             box(
                               title = "Venn Diagram", width = 8, status = "primary", #width = 6
                               plotOutput(outputId = "venn", height = "800px")
                             ),
                             box(
                               title = "Proteins identified", width = 6, status = "primary",
                               plotOutput(outputId = "protident", height = "600px")
                             )
                             
                           )
                         )),
                
                #Quality metrics tab
                
                tabPanel(title = "Quality control",
                         value = "quality",
                         icon = icon("chart-line"),
                         fluidPage(
                           fluidRow(
                             box(
                               title = "Dispersion plots", width = 6, status = "primary",
                               plotOutput(outputId = "boxplot", height = "600px")
                             ),
                             box(
                               title = "Imputation control plot", width = 6, status = "primary",
                               plotOutput(outputId = "preimputationplot", height = "600px")
                             ),
                             box(
                               title = "Normality plots", width = 6, status = "primary",
                               plotOutput(outputId = "histogram", height = "600px")
                             ),
                             box(
                               title = "Dimension reduction analysis", width = 6, status = "primary",
                               plotOutput(outputId = "pcaplot", height = "600px")
                             ),
                             box(
                               title = "Correlation plots", width = 8, status = "primary",
                               plotOutput(outputId = "sscatplot", height = "600px")
                             )
                             
                           )
                         )),
                
                #Differential analysis tab
                tabPanel(title = "Differential analysis",
                         value = "statistics",
                         icon = icon("table"),
                         fluidPage(
                           fluidRow(
                             infoBoxOutput("significantBox_dm",width = 4),
                             column( width = 12, DT::dataTableOutput("limma"))
                           )
                         )),
                
                #Differential analysis plots tab
                
                tabPanel(title = "Differential analysis plots",
                         value = "diffplots",
                         icon = icon("chart-area"),
                         fluidPage(
                           fluidRow(
                             box(
                               title = "Volcano plot", width = 6, status = "primary",
                               plotlyOutput(outputId = "volcanoplot", height = "600px")
                             ),
                             box(
                               title = "Heatmap", width = 6, status = "primary",
                               plotOutput(outputId = "heatmapplot", height = "600px")
                             )
                           ),
                           fluidRow(
                               box(
                               title = "Differential heatmap", width = 6, status = "primary",
                               plotOutput(outputId = "difheatmapplot", height = "600px")
                             ),
                             
                             box(
                               title = "Protein Intensity", width = 6, status = "primary",
                               plotOutput(outputId = "difboxplotplot", height = "600px")
                             ))
                             
                           
                           
                         ),
                ),
                
                #Functional analysis
                
                tabPanel( title = "Functional Analysis",
                          value = "fanalisis",
                          icon = icon("bezier-curve"),
                          column( width = 12, DT::dataTableOutput("functionalanalysis")),
                          fluidPage(
                            fluidRow(
                              box(
                                title = "Dotplot", width = 6, status = "primary",
                                plotOutput(outputId = "dotplotout", height = "600px")
                              ),
                              box(
                                title = "Barplot", width = 6, status = "primary",
                                plotOutput(outputId = "barplotout", height = "600px")
                              ),
                              box(
                                title = "Manhattan plot", width = 6, status = "primary",
                                plotlyOutput(outputId = "manhattanout", height = "600px")
                              ),
                            )
                          )),
                
                #Interactions
                
                tabPanel(title = "Interaction analysis",
                         value = "interactions",
                         icon = icon("atom"),
                         fluidPage(
                           fluidRow(
                             box(
                               title = "Up-regulated", width = 6, status = "primary",
                               plotOutput(outputId = "upreg", height = "600px")
                             ),
                             box(
                               title = "Down-regulated", width = 6, status = "primary",
                               plotOutput(outputId = "downreg", height = "600px")
                             ),
                             column(width = 6, DT::dataTableOutput("upgraph")),
                             column(width = 6, DT::dataTableOutput("downgraph")))
                           
                         )
                )
    )
  ))




server <- function(input, output) {
  #Home
  output$home_img <- renderImage({
    
    list(src = "www/traianprot.png",
         width = "80%",
         height = 1000)
    
  }, deleteFile = F)
  
  #Data handling
  #Several reactive variables in order to plot the quality metrics plots
  
  data_quick <- reactive({
    if (input$comptplatform == 1 | input$comptplatform == 2 | input$comptplatform == 5){
      req(input$file) #Controls whether a file has been uploaded or not
  
      raw <- read.delim(input$file$datapath, sep = "\t", stringsAsFactors = FALSE, colClasses = "character") 
      df<- quick_filtering(raw, input$intensitymode, input$comptplatform, input$organismfocus, input$condition1, input$condition2)
      return(df)
    } else if (input$comptplatform == 3){
      req(input$file) #Controls whether a file has been uploaded or not
      
      raw <- read.delim(input$file$datapath, sep = "\t", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE) 
      df<- quick_filtering(raw, input$intensitymode, input$comptplatform, input$organismfocus, input$condition1, input$condition2) #Marcarlo como dataframe 
      df <- as.data.frame(df)
      return(df)
    } else if (input$comptplatform == 4){
      
      raw <- as.data.frame(readxl::read_xlsx(input$file$datapath))
      
      df<- quick_filtering(raw, input$intensitymode, input$comptplatform, input$organismfocus, input$condition1, input$condition2)
      return(df)
    }
    
  })
  
  LOG2.names <- reactive({ #Set of samples in our dataset
    LOG2.names <- obtain_LOG.names(data_quick())
    return(LOG2.names)
  })
  
  cond.names <- reactive({ #subset of samples of interest in our dataset
    condition1_names <- grep(input$condition1, LOG2.names(), value = TRUE)
    condition2_names <- grep(input$condition2, LOG2.names(), value = TRUE)
    intensity_names <- c(condition1_names,condition2_names)
    return(intensity_names)
  })
  
  unique <- reactive({
    
    conditions <- c(input$condition1, input$condition2)
    replicas_condicion1 <- input$repcond1
    replicas_condicion2 <- input$repcond2
    if (input$labeltype == 2){ #We add this detail in order to get rid of the proteins that have NAs in all the measurements with TMT label. Only quantified proteins remain
     data_quick <-  data_quick()[rowSums(is.na(data_quick()[, cond.names()])) != length(cond.names()), ]
    }
    unique_proteins <- obtain_unique_proteins(data_quick(), conditions, LOG2.names(), replicas_condicion1, replicas_condicion2)
    if (input$uniquefilter == TRUE){
      
      unique_proteins <- lapply(unique_proteins, unique_peptides_filter,input$condition1, input$condition2, input$numberuniquepep, input$comptplatform)
      
    }
    
    return(unique_proteins)
  })
  
  unique_control <- reactive({
    unique_proteins <- as.data.frame(unique()[1], check.names = FALSE)
    return(unique_proteins)
  })
  
  unique_treatment <- reactive({
    unique_proteins <- as.data.frame(unique()[2], check.names = FALSE)
    return(unique_proteins)
  })
  
  data_filtered <- reactive({
    
    conditions <- c(input$condition1, input$condition2)
    min_count <- c(input$minfiltro, input$minfiltro)
    df.F <- filter_valids(data_quick(), unique(), conditions, min_count, at_least_one <- FALSE, LOG2.names(), input$labeltype)
    if (input$uniquefilter == TRUE){
      
      df.F <- unique_peptides_filter(df.F, input$condition1, input$condition2, input$numberuniquepep, input$comptplatform)
    }
    return(df.F)
    
  })
  
  data <- reactive({
    
    if (input$normchoice == TRUE){
      if (input$normoptions == "Median centering"){
        df.F <- median_centering(data_filtered(), LOG2.names())
      } else if (input$normoptions != "Median centering"){
          df.F <- normalization_func(data_filtered(), LOG2.names(), input$normoptions)
        }
    } else if (input$normchoice == FALSE){
      df.F <- data_filtered()
    }
    validate(need(!is.null(input$imputation),
                  "Please select an option"))
    if (input$imputation == 1) { 
      df.FNI <- impute_KNN_data(as.data.frame(df.F), LOG2.names(), k = 5)
    } else if (input$imputation == 2){
      df.FNI <- impute_data(as.data.frame(df.F), LOG2.names())
    } else if (input$imputation == 3){
      df.FNI <- df.F
    }
    return(df.FNI)
    
  })
  
  total_dataset <- reactive({
    total <- bind_rows(data(), unique_control())
    total_total <- bind_rows(total, unique_treatment())
    return(total_total)
  })
    
  

  output$file <- DT::renderDataTable({
    if (input$displaytable == 0){
      return(NULL)
    } else if (input$displaytable > 0){
     validate(need(!is.null(input$displayproteins),
                    "Please select an option"))
     if (input$displayproteins == 1){
  
    datatable(data(), options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15, 20), 
                                     scrollX = T,
                                     autoWidth = TRUE ))
     } else if (input$displayproteins == 2){
    
    datatable(unique_control(), options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15, 20), 
                                     scrollX = T,
                                     autoWidth = TRUE ))
    
     } else if (input$displayproteins == 3){
       
       datatable(unique_treatment(), options = list(pageLength = 5,
                                        lengthMenu = c(5, 10, 15, 20), 
                                        scrollX = T,
                                        autoWidth = TRUE ))
     }
    
    }})
  
  output$downloaddfcomm <- downloadHandler(
    
    filename = function() {
      paste(input$filenamedownloadcomm, Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.csv(total_dataset(), file)
    }
  )
  
  output$downloaddfcontroluniques <- downloadHandler(
    
    filename = function() {
      paste(input$filenamedownloadcontroluniques, Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.csv(unique_control(), file)
    }
  )
  
  output$downloaddftreatmentuniques <- downloadHandler(
    
    filename = function() {
      paste(input$filenamedownloadtreatmentuniques, Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      write.csv(unique_treatment(), file)
    }
  )
  
  
  output$LOG2.names <- renderTable({
    LOG2.names <- obtain_LOG.names(data_quick())
   
  }, striped = TRUE, align = "c", bordered = TRUE)
  
  

  #Proteins identified
  proteins_identified <- reactive({
    identify_proteins(data_quick(), cond.names(), input$comptplatform, input$repcond1, input$repcond2)
    
  })
  
  output$protident <- renderPlot({
  try(proteins_identified(), silent = TRUE)
  })
  
  #Venn plot
  
  venn_diagram_plot <- reactive({
    grid.draw(venn_diagram(data_filtered(), unique(), input$condition1venn, input$condition2venn, input$color1, input$color2))
    
  })
  output$venn <- renderPlot({
  try(venn_diagram_plot(), silent = TRUE)
    })
  
  output$downloadvenn <- downloadHandler(
    filename = function() {
      paste(input$venn_name, ".", input$vennextension, sep = "")
    },
    content = function(file) {
      
      tiff(file, width = 12, height = 10, units = "in", res = 400)
      grid.draw(venn_diagram(data_filtered(), unique(), input$condition1venn, input$condition2venn, input$color1, input$color2))
      dev.off()
    }
  )
  
  #Quality metrics 
  
  #Distribution plots
  boxplot_distribution <- reactive({
    boxplot_function(data(), cond.names(), cex.axis = 0.5)
    
  })
  
  output$boxplot <- renderPlot({
    if (input$displaydistplots == 1) {
      boxplot <- try(boxplot_distribution(), silent = TRUE)
    } else if (input$displaydistplots == 2) {
      try(plotCV2(data()[,cond.names()],  trend = TRUE, main = "Dispersion check", cex = 0.2, pch = 16, xlab="Average log-intensity", ylab=expression("Relative standard deviation")), silent = TRUE)
    }
    
  })
  
  #Imputation plots 
  pre_imp_plot <- reactive({
    preimputation_state(data_filtered(), cond.names())
  })
  
  output$preimputationplot <- renderPlot({
      if (input$displayimpplots == 1) {
        try(pre_imp_plot(), silent = TRUE)
      
    } else if (input$displayimpplots == 2) {
      try(postimputation_state(data(), cond.names()), silent = TRUE)
    }
   
  })
  

  #PCA dimreduction
#  PC_Analysis <- reactive({
#    my_pca <- pca(data(),cond.names(), input$repcond1, input$repcond2)
#  })
  
#  tsne_Analysis <- reactive({
#    my_tsne <- tsne(data(), cond.names(), input$repcond1, input$repcond2)
#  })
  
  output$pcaplot <- renderPlot({
    if (input$dimreduction == 1) {
      try(pca(data(),cond.names(), input$repcond1, input$repcond2), silent = TRUE)
    } else if (input$dimreduction == 2) {
      print("Running t-SNE...")  # Debugging message
      tsne(data(), cond.names(), input$repcond1, input$repcond2,  perplexity_num = input$perplexitytsne)
    }
  })
  
  
  output$boxplot <- renderPlot({
    if (input$displaydistplots == 1) {
      boxplot <- try(boxplot_distribution(), silent = TRUE)
    } else if (input$displaydistplots == 2) {
      try(plotCV2(data()[,cond.names()],  trend = TRUE, main = "Dispersion check", cex = 0.2, pch = 16, xlab="Average log-intensity", ylab=expression("Relative standard deviation")), silent = TRUE)
    }
    
  })
  
  #Scatter plot
  
  correlation_plot <- reactive({
    corrplot_function(data()[cond.names()], input$dispmethod)
  
  })
  
  output$histogram <- renderPlot({
    if (input$displaynormplots == 1) {
      histogram(data(), input$histsample, input$histcol, input$histtitle)
    } else if (input$displaynormplots == 2) {
      qqplot_function(data(), input$scatsample1, input$scatsample2, input$qqcolor)
    }
    
  })
    
  #Correlation plots

  output$sscatplot <- renderPlot({
    if (input$displaycorrplots == 1) {
      scatterplot_function(data(), input$scatsample1, input$scatsample2)
    } else if (input$displaycorrplots == 2) {
      correlation_plot()
    }
    
  })

  output$downloadqmplot <- downloadHandler(
    filename = function() {
      paste(input$qualitymetricsfilename, ".tiff" , sep = "") 
    },
    content = function(file) {
      if (input$qualitymetricsplotchoice == 1) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        boxplot_function(data(), cond.names(), cex.axis = 0.5)
        dev.off()
      } else if (input$qualitymetricsplotchoice == 2) {  
        tiff(file, width = 5, height = 4, units = "in", res = 300)
        plotCV2(
          data()[, cond.names()],
          trend = TRUE,
          main = "Dispersion check",
          cex = 0.2,
          pch = 16,
          xlab = "Average log-intensity",
          ylab = expression("Relative standard deviation"))
        dev.off()
      
      } else if (input$qualitymetricsplotchoice == 3) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        preimputation_state(data_filtered(), cond.names())
        dev.off()
      } else if (input$qualitymetricsplotchoice == 4) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        postimputation_state(data(), cond.names())
        dev.off()
      } else if (input$qualitymetricsplotchoice == 5) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        histogram(data(), input$histsample, input$histcol, input$histtitle)
        dev.off()
      } else if (input$qualitymetricsplotchoice == 6) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        qqplot_function(data(), input$scatsample1, input$scatsample2, input$qqcolor)
        dev.off()
      } else if (input$qualitymetricsplotchoice == 7) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        scatterplot_function(data(), input$scatsample1, input$scatsample2)
        dev.off()
      } else if (input$qualitymetricsplotchoice == 8) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        correlation_plot()
        dev.off()
      } else if (input$qualitymetricsplotchoice == 9) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        pca(data(),cond.names(), input$repcond1, input$repcond2)
        dev.off()
      }else if (input$qualitymetricsplotchoice == 10) {
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        print(tsne(data(),cond.names(), input$repcond1, input$repcond2, perplexity_num = input$perplexitytsne))
        dev.off()
      }
    }
  )
  
  
  #Differential expression
  
  difexpression <- reactive({
    validate(need(!is.null(input$pvaladj),
                  "Please select an option"))
    statistic_dataframe <- statistical_analysis(data(), LOG2.names(), input$displaytest, input$testid, input$repcond1, input$repcond2, input$condition1, input$condition2, input$LogFCup, input$LogFCdown, input$sigcutoff, input$pvaladj, input$statselected, unique(), input$proteins, input$PSMaware, input$comptplatform)
    
    statistic_dataframe <- statistic_dataframe[, !(names(statistic_dataframe) %in% "Protein_description")] #In order to avoid losing "Protein_description" to "Protein_description_X or Y " and turn it into non detectable column 
    merged <- merge(total_dataset(), statistic_dataframe, by = "Protein")
    
    if (input$PSMaware == TRUE){
      merged <- merged %>%
        select("Protein", "Protein_description", "logFC", "sca.P.Value", "sca.adj.pval", "expression", everything())
    } else if (input$PSMaware == FALSE){
      merged <- merged %>%
        select("Protein", "Protein_description", "logFC", "p.value", "adj.P.Val", "expression", everything())
    }
   
    
    merged <- merged[order(merged$expression),]
    row.names(merged) <- merged$Protein
    
    
    return(merged)
    
  })
  
  result.all <- reactive({
    result <- merge(data(), difexpression(), by = "Protein")
    return(result)
  })
  
 
  output$significantBox_dm <- renderInfoBox({
    num_total <- difexpression() %>%
      nrow()
    num_signif <- nrow(subset(difexpression(), difexpression()$expression != "Unchanged"))
    
    frac <- num_signif / num_total
    
    info_box <- 		infoBox("Significant proteins",
                          paste0(num_signif,
                                 " out of ",
                                 num_total),
                          paste0(signif(frac * 100, digits = 3),
                                 "% of proteins differentially expressed across all conditions"),
                          icon = icon("stats", lib = "glyphicon"),
                          color = "blue",
                          # fill = TRUE,
                          width = 10)
    
    return(info_box)
  })
  
  output$limma <- DT::renderDataTable({
      
      datatable(difexpression(), options = list(pageLength = 15,
                                                lengthMenu = c(5, 10, 15, 20), 
                                                scrollX = T,
                                                autoWidth = TRUE ))
      
    
  })
  
  output$downloaddif <- downloadHandler(
    filename = function() {
      paste(input$namedownload, Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(difexpression(), file)
    }
    
  )
  
  #Differential expression plots
  
  #Volcano
  volcano <- reactive({
    volcano <- volcano_plot(difexpression(), input$volcanotitle, input$labelprots, input$statselected, input$PSMaware)
    return(volcano)
  })
  
  volcano_tiff <- reactive({
    volcano <- volcano_plot_tiff(difexpression(), input$volcanotitle, input$statselected, input$PSMaware)
    return(volcano)
  })
  
  output$volcanoplot <- renderPlotly({
    volcano()
  })
  
  output$downloadvolcano <- downloadHandler(
    filename = function() {
      paste(input$name_download_volcano, ".", input$difextension, sep = "")
    },
    content = function(file) {
      if (input$difextension == "tiff"){
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        p <- volcano_plot_tiff(difexpression(), input$volcanotitle, input$statselected, input$PSMaware)
        print(p)
        dev.off()
      } else if (input$difextension == "pdf") {
        pdf(file, width = 12, height = 10, units = "in", res = 400)
        p <- volcano_plot_tiff(difexpression(), input$volcanotitle, input$statselected, input$PSMaware)
        print(p)
        dev.off()
      } else if (input$difextension == "jpeg"){
        jpeg(file, width = 12, height = 10, units = "in", res = 400)
        p <- volcano_plot_tiff(difexpression(), input$volcanotitle, input$statselected, input$PSMaware)
        print(p)
        dev.off()
      } else if (input$difextension == "png"){
        png(file, width = 12, height = 10, units = "in", res = 400)
        p <- volcano_plot_tiff(difexpression(), input$volcanotitle, input$statselected, input$PSMaware)
        print(p)
        dev.off()
      }
    }
  )
  
  
  #Heatmap
  heatmap <- reactive({
    heatmap <- my_heatmap(data(), cond.names(), input$heatmaptitle)
    return(heatmap)
  })
  
  output$heatmapplot <- renderPlot({
    tryCatch({
      # Your normal plot code here
      heatmap(your_data)
      
    }, error = function(e) {
      # Show a blank plot with a custom message
      plot.new()
      text(0.5, 0.5, "Something went wrong while generating the plot, check whether the amount of NAs is too large.", cex = 1.2)
      # Optionally log the original error for debugging:
      # message("Plot error: ", e$message)
    })
  })
  
  output$downloadheatmap <- downloadHandler(
    filename = function() {
      paste(input$name_download_heatmap, ".", input$difextension, sep = "")
    },
    content = function(file) {
      if (input$difextension == "tiff"){
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap(data(), cond.names(), input$heatmaptitle)
        dev.off()
      } else if (input$difextension == "pdf") {
        pdf(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap(data(), cond.names(), input$heatmaptitle)
        dev.off()
      } else if (input$difextension == "jpeg"){
        jpeg(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap(data(), cond.names(), input$heatmaptitle)
        dev.off()
      } else if (input$difextension == "png"){
        png(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap(data(), cond.names(), input$heatmaptitle)
        dev.off()
      }
    }
  )
  
 
  
  #Differential heatmap
  
  diff_heatmap <- reactive({
    diff_heatmap <- my_heatmap_differential(difexpression(),data(), cond.names(), input$heatmaptitle)
    return(diff_heatmap)
  })
  
  output$difheatmapplot <- renderPlot({
    tryCatch({
      # Your normal plot code here
      diff_heatmap()
      
    }, error = function(e) {
      # Show a blank plot with a custom message
      plot.new()
      text(0.5, 0.5, "Please check whether the amount of NAs is too large.", cex = 1.2)
      # Optionally log the original error for debugging:
      # message("Plot error: ", e$message)
    })
  })

  output$downloaddifheatmap <- downloadHandler(
    filename = function() {
      paste(input$name_download_heatmap, ".", input$difextension, sep = "")
    },
    content = function(file) {
      if (input$difextension == "tiff"){
        tiff(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap_differential(difexpression(),data(), cond.names(), input$heatmaptitle)
        dev.off()
      } else if (input$difextension == "pdf") {
        pdf(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap_differential(difexpression(),data(), cond.names(), input$heatmaptitle)
        dev.off()
      } else if (input$difextension == "jpeg"){
        jpeg(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap_differential(difexpression(),data(), cond.names(), input$heatmaptitle)
        dev.off()
      } else if (input$difextension == "png"){
        png(file, width = 12, height = 10, units = "in", res = 400)
        my_heatmap_differential(difexpression(),data(), cond.names(), input$heatmaptitle)
        dev.off()
      }
    }
  )
  
  #Differential Boxplot
  dif_boxplot <- reactive({
    
    dif_boxplot <- Diferential_boxplot(data(), first_condition = input$condition1, second_condition = input$condition2, protein = input$proteinname, LOG2.names())
    
    return(dif_boxplot)
    
  })
  
  output$difboxplotplot <- renderPlot({
    dif_boxplot()
  })
  
  
  #Functional analysis 
  
  funcanalysis <- reactive({
    validate(need(!is.null(input$organism),
                  "Please a valid organism"))
    validate(need(!is.null(input$proteinid),
                  "Please a valid Protein ID"))
    Go_terms <- Goterms_finder(difexpression(), data_quick(), target = input$proteinid, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organismo = input$organism, custombg = input$backgroundset, input$comptplatform,  user_threshold = input$userthreshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC"))
    return(Go_terms)  
  })
  
  output$functionalanalysis <- DT::renderDataTable({
    
    datatable(funcanalysis()[[2]]@result, options = list(pageLength = 10,
                                                         lengthMenu = c(5, 10, 15, 20), 
                                                         scrollX = T,
                                                         autoWidth = TRUE ))
  })
  
  output$downloadfunc <- downloadHandler(
    filename = function() {
      paste(input$name_download, Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      writexl::write_xlsx(funcanalysis()[[2]]@result, file)
      
    }
  )
  
  #Dotplot
  dotplot_react <- reactive({
    dotplot_func(funcanalysis(), x = "GeneRatio", title = input$dotplottitle, split = "Conditions", font.size = input$fontsizedotplot, showCategory = input$categorydotplot, color = "adj.P.Val")
    
  })
  
  output$dotplotout <- renderPlot({
    dotplot_react()
  })
  
  
  output$downloaddotplot <- downloadHandler(
    filename = function() {
      paste(input$name_download_func, ".", input$funcextension, sep = "")
    },
    content = function(file) {
      ggsave( file,
              plot = dotplot_func(funcanalysis(), x = "GeneRatio", title = input$dotplottitle, split = "Conditions", font.size = input$fontsizedotplot, showCategory = input$categorydotplot, color = "adj.P.Val"),
              device = input$funcextension,
              dpi = input$funcquality)
    }
  )
  
  #Barplot
  barplot_react <- reactive({
    barplot_func(funcanalysis(), input$categorybarplot, conditions = input$barplottitle, font.size = input$fontsizebarplot)
  })
  
  output$barplotout <- renderPlot({
    barplot_react()
  })
  
  
  output$downloadbarplot <- downloadHandler(
    filename = function() {
      paste(input$name_download_func, ".", input$funcextension, sep = "")
    },
    content = function(file) {
      ggsave( file,
              plot =  barplot_func(funcanalysis(), input$categorybarplot, conditions = input$barplottitle, font.size = input$fontsizebarplot),
              device = input$funcextension,
              dpi = input$funcquality)
    }
  )
  
  #Manhattan plot
  manhattan <- reactive({
    manhattan <- gostplot_func(funcanalysis())
    return(manhattan)
  })
  
  output$manhattanout <- renderPlotly({
    if (input$manhattan == 0){
      return(NULL)
    } else if (input$manhattan > 0){
      manhattan()
    }
  })
  
  #Interaction analysis
  
  interactions <- reactive({
    interactions <- interactions_up(difexpression(), input$taxonid, input$scthreshold)
    return(interactions)
    
  })
  

  
  output$upreg <- renderPlot({
    if (input$upnet == 0){
      return(NULL)
    } else if (input$upnet > 0){
      interactions_up(difexpression(), input$taxonid, input$scthreshold)
    }
  }
  )
  
  output$downloadupnetwork <- downloadHandler(
    filename = function() {
      paste(input$name_download_int, ".", input$intextension, sep = "")
    },
    content = function(file) {
      if (input$funcextension == "tiff"){
        tiff(file)
        interactions_up(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      } else if (input$funcextension == "pdf") {
        pdf(file)
        interactions_up(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      } else if (input$funcextension == "jpeg"){
        jpeg(file)
        interactions_up(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      } else if (input$funcextension == "png"){
        png(file)
        interactions_up(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      }
    }
  )
  
  output$downreg <- renderPlot({
    if (input$downnet == 0){
      return(NULL)
    } else if (input$downnet > 0){
      interactions_down(difexpression(), input$taxonid, input$scthreshold)
    }
  }
  
  )
  
  output$downloaddownnetwork <- downloadHandler(
    filename = function() {
      paste(input$name_download_int, ".", input$intextension, sep = "")
    },
    content = function(file) {
      if (input$funcextension == "tiff"){
        tiff(file)
        interactions_down(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      } else if (input$funcextension == "pdf") {
        pdf(file)
        interactions_down(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      } else if (input$funcextension == "jpeg"){
        jpeg(file)
        interactions_down(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      } else if (input$funcextension == "png"){
        png(file)
        interactions_down(difexpression(), input$taxonid, input$scthreshold)
        dev.off()
      }
    }
  )
  
  #Graph analysis 
  
  graph <- reactive({
    graph_analysis <- igraph_analysis(interactions(), input$taxonid, input$scthreshold)
    return(graph_analysis)
  })
  
  output$upgraph <- DT::renderDataTable({

    datatable(graph()[[1]], options = list(pageLength = 1,
                                     lengthMenu = c(1), 
                                     scrollX = T,
                                     autoWidth = TRUE ))
    
  })
  
  output$downgraph <- DT::renderDataTable({
    
    datatable(graph()[[2]], options = list(pageLength = 1,
                                           lengthMenu = c(1), 
                                           scrollX = T,
                                           autoWidth = TRUE ))
    
  })
  
  ##### ====== Download report ====== #####
  
  output$downloadReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Downstreaming_Analysis_report.Rmd")
      file.copy("../Downstreaming_Analysis_report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        venn = venn_diagram_plot(),
        proteins = proteins_identified(),
        boxplot = boxplot_distribution(),
        imputation = pre_imp_plot(),
        pca = PC_Analysis(),
        corrplot = correlation_plot(),
        volcano = volcano(),
        heatmap = heatmap(),
        dotplot = dotplot_react(),
        barplot = barplot_react(),
        control_replicates = input$repcond1,
        treatment_replicates = input$repcond2,
        type_of_test = input$testid,
        log2FC_up = input$LogFCup,
        log2FC_down = input$LogFCdown,
        pvalue = input$sigcutoff,
        organism = input$organism,
        p_value_for_enrichment = input$userthreshold,
        taxon_id = input$taxonid,
        score_threshold = input$scthreshold,
      )
      
      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport, output_format = "html_document", output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  

}

shinyApp(ui, server)

