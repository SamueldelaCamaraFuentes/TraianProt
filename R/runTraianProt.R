#'
#' @title Launch the TraianProt Shiny Application
#' @description This function builds and returns the Shiny app object for TraianProt.
#'
#' @return A \code{shiny.appobj} object representing the TraianProt app.
#' @export
#' @import shiny
#' @import shinydashboard
#' @importFrom shinyWidgets switchInput
#' @importFrom DT renderDT DTOutput
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom dplyr mutate select filter pull bind_rows %>%
#'
#' @examples
#' app <- runTraianProt()
#'
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
runTraianProt <- function() {
    ui <- .build_ui()
    server <- .build_server

    app <- shiny::shinyApp(ui = ui, server = server)
    return(app)
}


# ==========================================
# UI BUILDER & SUB-COMPONENTS
# ==========================================

.build_ui <- function() {
    dashboardPage(
        skin = "red",
        dashboardHeader(title = "TraianProt"),
        sidebar = .build_sidebar(),
        body = .build_body()
    )
}


.build_sidebar <- function() {
    dashboardSidebar(
        # Home Panel
        conditionalPanel(
            condition = "input.main_tabs == 'Home_introduction'",
            sidebarMenu(
                menuItem("Home",
                    icon = icon("home"), selected = TRUE, tabName = "home",
                    downloadButton("download_tutorial", "Download User Guide",
                        icon = icon("file-pdf"),
                        style = "display: block; margin: 10px auto; width: 200px; color:black;"
                    )
                )
            )
        ),

        # Data Handling Panel
        conditionalPanel(
            condition = "input.main_tabs == 'data_handling'",
            sidebarMenu(
                menuItem("File Input",
                    tabName = "file input", icon = icon("upload"),
                    fileInput(inputId = "file", label = h3("File input"), multiple = FALSE, accept = c(".csv", ".tsv", ".txt")),
                    fileInput(inputId = "metadata", label = h3("Metadata"), multiple = FALSE, accept = c(".csv", ".tsv", ".txt")),
                    switchInput(inputId = "technicalreplicates", label = "Technical replicates", value = FALSE, onLabel = "Yes", offLabel = "No"),
                    selectInput(inputId = "comptplatform", label = "Platform", choices = list("MaxQuant" = 1, "MSFragger" = 2, "DIA-NN" = 3, "Proteome Discoverer" = 4, "ProteoScape" = 5), selected = 1),
                    conditionalPanel(condition = "input.comptplatform == 3", textInput(inputId = "pathdir", label = "Directory Path:", value = getwd())),
                    actionButton("setdir", "Set Directory"),
                    selectInput(inputId = "labeltype", label = "Type", choices = list("Label free" = 1, "TMT" = 2, "SILAC" = 3), selected = 1),
                    uiOutput("pairwise_selector"),
                    selectInput(inputId = "organismfocus", label = "Organism", choices = list("Candida albicans" = 1, "Other" = 2), selected = 1)
                ),
                menuItem("Pre-processing",
                    icon = icon("th"), tabName = "preprocessing",
                    h4("Filtering by samples"), numericInput("min_prop", "Proportion for filtering", value = 0.5),
                    h4("Additional filtering step"), numericInput("numberuniquepep", "Minimum number of unique peptides", value = 1),
                    numericInput("proportionsamples", "Proportion of samples", value = 0.5),
                    switchInput("uniquefilter", "Unique peptides filter", value = TRUE, onLabel = "Yes", offLabel = "No"),
                    h3("Normalization"), switchInput("normchoice", "Normalization", value = TRUE, onLabel = "Yes", offLabel = "No"),
                    textInput("normoptions", "Normalization method", value = "Median centering"), helpText("Note: mean, median, trimMean and vsc"),
                    h3("Imputation"), selectInput("imputation", "Imputation method", choices = list("No imputation" = 3, "K Nearest Neighbors" = 1, "Normal distribution" = 2), selected = 3),
                    actionButton("displaytable", "Display table"),
                    selectInput("displayproteins", "Choose a protein set to display", choices = list("Common proteins between conditions" = 1, "Exclusive control proteins" = 2, "Exclusive treatment proteins" = 3), selected = 1),
                    downloadButton("downloadprotquant", "Proteins Quantified", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;")
                ),
                menuItem("Venn Diagram",
                    icon = icon("th"), tabName = "venn",
                    strong(h4("Condition 1")), textInput("color1", "Color:", value = "blue"),
                    strong(h4("Condition 2")), textInput("color2", "Color:", value = "maroon"),
                    selectInput("vennextension", "File type:", choices = c("tiff", "jpeg", "png"), selected = "pdf"),
                    selectInput("vennquality", "Quality:", choices = c("High" = "retina", "Medium" = "print", "low" = "screen"), selected = ""),
                    textInput("venn_name", "Filename Venn", value = "Venn"), h5("Venn Diagram download"),
                    downloadButton("downloadvenn", "Download Venn", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;")
                ),
                menuItem("Power of a test",
                    icon = icon("th"), tabName = "powertest",
                    strong(h4("Fold Change")), numericInput("FoldChangepower", "Fold change", value = 2),
                    numericInput("replicatespower", "Maximun number of replicates", value = 30),
                    selectInput("powerchoice", "Choice:", choices = list("p value" = 1, "adjusted p value" = 2, selected = 2)),
                    numericInput("statisticalpower", "value", value = 0.05)
                ),
                menuItem("Unique Peptides Extractor", icon = icon("external-link-alt"), href = "https://samueldelacamara.shinyapps.io/Unique_peptides_extractor/", newtab = TRUE),
                menuItem("Download tables",
                    icon = icon("download"), tabName = "Download",
                    textInput("filenamedownloadcomm", "Filename proteins", value = "data"), downloadButton("downloaddfcomm", icon("download"), label = "Download", style = "display: block; margin: 0 auto; width: 200px; color:black;"),
                    textInput("filenamedownloadcontroluniques", "Filename unique control proteins", value = "data"), downloadButton("downloaddfcontroluniques", icon("download"), label = "Download", style = "display: block; margin: 0 auto; width: 200px; color:black;"),
                    textInput("filenamedownloadtreatmentuniques", "Filename treatment unique proteins", value = "data"), downloadButton("downloaddftreatmentuniques", icon("download"), label = "Download", style = "display: block; margin: 0 auto; width: 200px; color:black;")
                )
            )
        ),

        # Quality Panel Sidebar
        conditionalPanel(
            condition = "input.main_tabs == 'quality'",
            sidebarMenu(
                menuItem("Quality metrics",
                    icon = icon("chart-line"),
                    strong(h4("Distribution plots")), selectInput("displaydistplots", "Select a plot", choices = list("Boxplot" = 1, "Dispersion plot" = 2), selected = 1),
                    strong(h4("Imputation plot")), selectInput("displayimpplots", "Select a plot", choices = list("Pre imputation" = 1, "Post imputation" = 2), selected = 1),
                    strong(h4("Normality plots")), textInput("histsample", "Sample for histogram:", value = "LOG2.WT1"), textInput("histcol", "Color for histogram:", value = "maroon"), textInput("histtitle", "Title for histogram:", value = "Histogram of intensities"),
                    textInput("qqcolor", "Color for QQplot:", value = "blue"), selectInput("displaynormplots", "Select a plot", choices = list("Histogram" = 1, "Q-Q plot" = 2), selected = 1),
                    strong(h4("Dimension reduction")), selectInput("dimreduction", "Select a plot", choices = list("PCA" = 1, "t-SNE" = 2), selected = 1),
                    conditionalPanel(condition = "input.dimreduction == '1'", selectInput("pc_x_axis", "Select PC for X-axis:", choices = seq_len(5), selected = 1), selectInput("pc_y_axis", "Select PC for Y-axis:", choices = seq_len(5), selected = 2)),
                    numericInput("perplexitytsne", "perplexity tsne", value = 1),
                    strong(h4("Correlation plots")), textInput("scatsample1", "Sample1 for scatter:", value = "LOG2.WT1"), textInput("scatsample2", "Sample2 for scatter:", value = "LOG2.WT2"), textInput("dispmethod", "Display method:", value = "shade"),
                    selectInput("displaycorrplots", "Select a plot", choices = list("Scatter plot" = 1, "Correlation plot" = 2), selected = 1)
                ),
                menuItem("Download plots",
                    icon = icon("download"), tabName = "Download", strong(h4("Quality metrics download")),
                    selectInput("qualitymetricsplotchoice", "Select a plot", choices = list("Boxplot" = 1, "Dispersion plot" = 2, "Pre imputation" = 3, "Post imputation" = 4, "Histogram" = 5, "Q-Q plot" = 6, "Scatter plot" = 7, "Correlation plot" = 8, "PCA" = 9, "t-SNE" = 10), selected = 1),
                    textInput("qualitymetricsfilename", "Filename", value = "quality metrics plot"), downloadButton("downloadqmplot", icon("download"), label = "Download", style = "display: block; margin: 0 auto; width: 200px; color:black;")
                )
            )
        ),

        # Statistics Panel Sidebar
        conditionalPanel(
            condition = "input.main_tabs == 'statistics'",
            sidebarMenu(
                menuItem("Differential analysis",
                    icon = icon("table"),
                    h3("Tipe of test"), selectInput("displaytest", "Choose a test to analyze:", choices = list("Simple t test approach" = 1, "Limma approach" = 2, "Wilcoxon test" = 3), selected = 2),
                    switchInput("PSMaware", "PSMs correction", value = TRUE, onLabel = "Yes", offLabel = "No"),
                    selectInput("proteins", "Choose a way to proceed:", choices = list("Common proteins" = 1, "All protein" = 2), selected = 2),
                    switchInput("testid", "Test", value = FALSE, onLabel = "Paired t-test", offLabel = "Two sample t-test"),
                    h3("Choice:"), selectInput("statselected", "Select a statistic:", choices = list("P value" = 1, "Q value" = 2), selected = 2),
                    numericInput("LogFC", "Log2FC threshold:", value = 1), numericInput("sigcutoff", "Statistic threshold:", value = 0.05),
                    selectInput("pvaladj", "Choose a p-value adjustment:", choices = list("FDR" = "fdr", "Bonferroni" = "bonferroni", "BH" = "BH", "Hochberg" = "hochberg"), selected = "fdr")
                ),
                menuItem("Download tables", icon = icon("download"), tabName = "Download", textInput("namedownload", "Filename", value = "differential data"), downloadButton("downloaddif", icon("download"), label = "Download", style = "display: block; margin: 0 auto; width: 200px; color:black;"))
            )
        ),

        # Differential Plots Panel Sidebar
        conditionalPanel(
            condition = "input.main_tabs == 'diffplots'",
            sidebarMenu(
                menuItem("Differential analysis plots", icon = icon("chart-area"), strong(h4("Volcano plot")), numericInput("labelprots", "Point size", value = 10), textInput("volcanotitle", "Insert a title", value = "Treatment vs Control"), strong(h4("Heatmap")), textInput("heatmaptitle", "Insert a title", value = "Treatment vs Control"), strong(h4("Protein Intensity")), textInput("proteinname", "Insert a title", value = "Treatment vs Control")),
                menuItem("Download plot options", icon = icon("download"), selectInput("difextension", "File type:", choices = c("tiff", "pdf", "jpeg", "png"), selected = "pdf"), selectInput("difquality", "Quality:", choices = c("High" = "retina", "Medium" = "print", "low" = "screen"), selected = ""), textInput("name_download_volcano", "Filename Volcano", value = "Volcano"), textInput("name_download_heatmap", "Filename Heatmap", value = "Heatmap"), h5("Volcano download"), downloadButton("downloadvolcano", "Download volcano", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;"), h5("Heatmap download"), downloadButton("downloadheatmap", "Download heatmap", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;"), h5("Differential Heatmap download"), downloadButton("downloaddifheatmap", "Download heatmap", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;"))
            )
        ),

        # Functional Analysis Panel Sidebar
        conditionalPanel(
            condition = "input.main_tabs == 'fanalisis'",
            sidebarMenu(
                menuItem("Functional analysis", icon = icon("table"), textInput("organism", "Organism:", value = "calbicans"), textInput("proteinid", "Type of Protein ID:", value = "ENSG"), numericInput("userthreshold", "enrichment threshold:", value = 0.05), switchInput("backgroundset", "Background", value = FALSE, onLabel = "Yes", offLabel = "No")),
                menuItem("Plots", icon = icon("chart-area"), strong(h4("Dotplot")), textInput("dotplottitle", "Title:", value = "Treatment vs Control"), numericInput("categorydotplot", "Number of terms to display:", value = 10), numericInput("fontsizedotplot", "Font size:", value = 16), strong(h4("Barplot")), textInput("barplottitle", "Title:", value = "Treatment vs Control"), numericInput("categorybarplot", "Number of terms to display:", value = 10), numericInput("fontsizebarplot", "Font size:", value = 10), strong(h4("Manhattan plot")), actionButton("manhattan", "Render plot")),
                menuItem("Download plot options", icon = icon("download"), selectInput("funcextension", "File type:", choices = c("tiff", "pdf", "jpeg", "png"), selected = "pdf"), selectInput("funcquality", "Quality:", choices = c("High" = "retina", "Medium" = "print", "low" = "screen"), selected = ""), textInput("name_download_func", "Filename", value = "functional analysis data"), h5("Dotplot download"), downloadButton("downloaddotplot", "Download", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;"), h5("Barplot download"), downloadButton("downloadbarplot", "Download", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;")),
                menuItem("Download tables", icon = icon("download"), tabName = "Download", textInput("name_download", "Filename", value = "functional analysis data"), downloadButton("downloadfunc", icon("download"), label = "Download", style = "display: block; margin: 0 auto; width: 200px; color:black;"))
            )
        ),

        # Interactions Panel Sidebar
        conditionalPanel(
            condition = "input.main_tabs == 'interactions'",
            sidebarMenu(
                menuItem("STRINGdb", icon = icon("atom"), strong(h4("Species taxon id")), numericInput("taxonid", "Taxon id", value = 237561), strong(h4("Score threshold")), numericInput("scthreshold", "Score threshold", value = 400), strong(h4("Up-regulated network")), actionButton("upnet", "Render plot"), strong(h4("Down-regulated network")), actionButton("downnet", "Render plot")),
                menuItem("Download plot options", icon = icon("download"), selectInput("intextension", "File type:", choices = c("tiff", "pdf", "jpeg", "png"), selected = "pdf"), selectInput("intquality", "Quality:", choices = c("High" = "retina", "Medium" = "print", "low" = "screen"), selected = ""), textInput("name_download_int", "Filename", value = "functional analysis data"), h5("Up-regulated download"), downloadButton("downloadupnetwork", "Download", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;"), h5("Down-regulated download"), downloadButton("downloaddownnetwork", "Download", icon = icon("download"), style = "display: block; margin: 0 auto; width: 200px; color:black;"))
            )
        )
    )
}

.build_body <- function() {
    dashboardBody(
        tabsetPanel(
            id = "main_tabs",
            .tab_home(),
            .tab_preprocessing(),
            .tab_quality_control(),
            .tab_differential_analysis(),
            .tab_differential_plots(),
            .tab_functional_analysis(),
            .tab_interaction_analysis(),
            .tab_report()
        )
    )
}

# Sub-tab UI functions
.tab_home <- function() {
    tabPanel(
        title = "Home", value = "Home_introduction", icon = icon("home"),
        div(imageOutput("home_img", inline = TRUE), style = "text-align: center; margin-left: 0%;"),
        div(
            style = "border-radius: 10px; background-color: #f9f9f9; border: 2px solid #ddd; padding: 20px; margin-top: 20px; width: 50%; box-shadow: 0px 4px 10px rgba(0,0,0,0.1); text-align: center; margin-left: auto; margin-right: auto; font-family: 'Arial', sans-serif;",
            h4("Welcome to TraianProt", style = "color: #2c3e50; font-weight: bold;"),
            p("TraianProt is a web-based, user-friendly proteomics data analysis platform...", style = "color: #34495e; font-size: 16px;")
        ),
        br(), hr(),
        h3("Try it out!", style = "text-align: center;"),
        p("Don't have data handy? Download our sample datasets.", style = "text-align: center;"),
        fluidRow(
            column(
                width = 12, offset = 3,
                box(
                    title = "Example Data for Testing", status = "danger", solidHeader = TRUE, width = 6,
                    "Download these files and upload them in the 'File Input' tab.", br(), br(),
                    downloadButton("dl_example_matrix", "Download ProteinGroups Example (MaxQuant)", style = "color: #fff; background-color: #d9534f; border-color: #d43f3a; width: 100%;"), br(), br(),
                    downloadButton("dl_example_metadata", "Download Metadata Example (MaxQuant)", style = "color: #fff; background-color: #d9534f; border-color: #d43f3a; width: 100%;")
                )
            )
        )
    )
}

.tab_preprocessing <- function() {
    tabPanel(
        title = "Preprocessing", value = "data_handling", icon = icon("table"),
        box(title = "Note", width = 12, status = "warning", solidHeader = TRUE, helpText("Tip: For FragPipe and Proteome Discoverer datasets...")),
        column(width = 12, DT::dataTableOutput("file")),
        fluidPage(
            fluidRow(
                column(3, h3("Columns"), tableOutput("LOG2.names")),
                box(title = "Venn Diagram", width = 8, status = "primary", plotOutput("venn", height = "800px")),
                box(title = "Proteins identified", width = 6, status = "primary", plotOutput("protident", height = "600px")),
                box(title = "Power of a test", width = 6, status = "primary", plotOutput("powertest", height = "600px")),
                box(title = "Note", width = 12, status = "warning", solidHeader = TRUE, helpText("If your data originates from FragPipe..."))
            )
        )
    )
}

.tab_quality_control <- function() {
    tabPanel(
        title = "Quality control", value = "quality", icon = icon("chart-line"),
        fluidPage(
            fluidRow(
                box(title = "Dispersion plots", width = 6, status = "primary", plotOutput("boxplot", height = "600px")),
                box(title = "Imputation control plot", width = 6, status = "primary", plotOutput("preimputationplot", height = "600px")),
                box(title = "Normality plots", width = 6, status = "primary", plotOutput("histogram", height = "600px")),
                box(title = "Dimension reduction analysis", width = 6, status = "primary", plotOutput("pcaplot", height = "600px")),
                box(title = "Correlation plots", width = 8, status = "primary", plotOutput("sscatplot", height = "600px"))
            )
        )
    )
}

.tab_differential_analysis <- function() {
    tabPanel(
        title = "Differential analysis", value = "statistics", icon = icon("table"),
        fluidPage(fluidRow(infoBoxOutput("significantBox_dm", width = 4), column(width = 12, DT::dataTableOutput("limma"))))
    )
}

.tab_differential_plots <- function() {
    tabPanel(
        title = "Differential analysis plots", value = "diffplots", icon = icon("chart-area"),
        fluidPage(
            fluidRow(
                box(title = "Volcano plot", width = 6, status = "primary", plotlyOutput("volcanoplot", height = "600px")),
                box(title = "Heatmap", width = 6, status = "primary", plotOutput("heatmapplot", height = "600px"))
            ),
            fluidRow(
                box(title = "Differential heatmap", width = 6, status = "primary", plotOutput("difheatmapplot", height = "600px")),
                box(title = "Protein Intensity", width = 6, status = "primary", plotOutput("difboxplotplot", height = "600px"))
            )
        )
    )
}

.tab_functional_analysis <- function() {
    tabPanel(
        title = "Functional Analysis", value = "fanalisis", icon = icon("bezier-curve"),
        column(width = 12, DT::dataTableOutput("functionalanalysis")),
        fluidPage(fluidRow(
            box(title = "Dotplot", width = 6, status = "primary", plotOutput("dotplotout", height = "600px")),
            box(title = "Barplot", width = 6, status = "primary", plotOutput("barplotout", height = "600px")),
            box(title = "Manhattan plot", width = 6, status = "primary", plotlyOutput("manhattanout", height = "600px"))
        ))
    )
}

.tab_interaction_analysis <- function() {
    tabPanel(
        title = "Interaction analysis", value = "interactions", icon = icon("atom"),
        fluidPage(fluidRow(
            box(title = "Up-regulated", width = 6, status = "primary", plotOutput("upreg", height = "600px")),
            box(title = "Down-regulated", width = 6, status = "primary", plotOutput("downreg", height = "600px")),
            column(width = 6, DT::dataTableOutput("upgraph")),
            column(width = 6, DT::dataTableOutput("downgraph"))
        ))
    )
}

.tab_report <- function() {
    tabPanel(
        title = "Report", icon = icon("file-alt"), value = "report_tab",
        div(
            style = "text-align: center; margin-top: 50px;",
            h1(icon("clipboard-check"), style = "color: #d9534f; font-size: 80px;"),
            h2("Analysis Complete!"), p("Generate a comprehensive HTML report...", style = "font-size: 18px; color: #555;"), br(),
            downloadButton("downloadreport", "Download Full Analysis Report", style = "color: #fff; background-color: #d9534f; border-color: #d43f3a; font-size: 20px; padding: 15px 30px;")
        )
    )
}


# ==========================================
# SERVER BUILDER & RESPONSIBILITY ROUTERS
# ==========================================

.build_server <- function(input, output, session) {
    # Max request size handling
    old_options <- options(shiny.maxRequestSize = 30 * 1024^2)
    session$onSessionEnded(function() {
        options(old_options)
    })

    # ----------------------------------------
    # CORE DATA PIPELINE (Shared Reactives)
    # ----------------------------------------

    metadata <- reactive({
        req(input$metadata$datapath)
        meta <- read.delim(input$metadata$datapath, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
        if (input$comptplatform %in% c(1, 2, 5)) {
            req(input$file)
            meta <- meta %>% mutate(raw_name = intensity_sample_name, intensity_sample_name = raw_name, log2_col = sub("Intensity", "LOG2", raw_name), unique_peptides_col = sub("Intensity", "Unique.peptides", raw_name))
        } else if (input$comptplatform == 3) {
            meta <- meta %>% mutate(raw_name = basename(intensity_sample_name), intensity_sample_name = raw_name, log2_col = sub("\\.(d|raw)$", ".LOG2", raw_name, ignore.case = TRUE), unique_peptides_col = paste0("Unique peptides ", sub("\\.(d|raw)$", "", raw_name, ignore.case = TRUE)))
        } else if (input$comptplatform == 4) {
            meta <- meta %>% mutate(raw_name = intensity_sample_name, intensity_sample_name = raw_name, log2_col = sub("Abundance:", "LOG2", raw_name), unique_peptides_col = sub("Abundance:", "Unique.peptides", raw_name))
        }
        return(meta %>% select(intensity_sample_name, group, sample_name, log2_col, unique_peptides_col))
    })
    condition_pairs <- reactive({
        req(metadata())
        conds <- unique(metadata()$group)
        if (length(conds) < 2) {
            return(NULL)
        }
        pairs <- expand.grid(conds, conds, stringsAsFactors = FALSE) %>% filter(Var1 != Var2)
        split(apply(pairs, 1, identity), rep(seq_len(nrow(pairs)), each = 2))
    })

    output$pairwise_selector <- renderUI({
        req(condition_pairs())
        # vapply devuelve character(1) y USE.NAMES = FALSE nos ahorra usar unname()
        selectInput("selected_pair", "Select Comparison", choices = vapply(condition_pairs(), function(p) paste(p[1], "vs", p[2]), character(1), USE.NAMES = FALSE))
    })

    selected_conditions <- reactive({
        req(condition_pairs())
        if (is.null(input$selected_pair)) unique(metadata()$group)[seq_len(2)] else strsplit(input$selected_pair, " vs ")[[1]]
    })

    filtered_metadata_base <- reactive({
        req(metadata(), selected_conditions())
        metadata() %>% filter(group %in% selected_conditions())
    })
    selected_dir <- reactiveVal(getwd())

    observeEvent(input$setdir, {
        req(input$pathdir)
        if (dir.exists(input$pathdir)) {
            selected_dir(normalizePath(input$pathdir))
        } else {
            showNotification("Error: Directory does not exist.", type = "error")
        }
    })

    data_quick_base <- reactive({
        req(input$file)
        if (input$comptplatform %in% c(1, 2, 5)) {
            raw <- read.delim(input$file$datapath, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
            quick_filtering(raw, input$comptplatform, input$organismfocus, filtered_metadata_base(), selected_conditions())
        } else if (input$comptplatform == 3) {
            raw <- read.delim(input$file$datapath, sep = "\t", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
            as.data.frame(quick_filtering(raw, input$comptplatform, input$organismfocus, filtered_metadata_base(), selected_conditions(), selected_dir()))
        } else if (input$comptplatform == 4) {
            raw <- as.data.frame(readxl::read_xlsx(input$file$datapath))
            quick_filtering(raw, input$comptplatform, input$organismfocus, filtered_metadata_base(), selected_conditions())
        }
    })

    collapsed_bundle <- reactive({
        req(data_quick_base(), filtered_metadata_base())
        collapse_technical_replicates(data_quick_base(), filtered_metadata_base())
    })
    data_quick <- reactive({
        if (input$technicalreplicates) collapsed_bundle()$data else data_quick_base()
    })
    filtered_metadata <- reactive({
        if (input$technicalreplicates) collapsed_bundle()$metadata else filtered_metadata_base()
    })

    data_quick_upf <- reactive({
        if (input$uniquefilter) unique_peptides_filter(data_quick(), filtered_metadata(), input$numberuniquepep, input$proportionsamples) else data_quick()
    })

    LOG2.names <- reactive({
        req(data_quick())
        obtain_LOG.names(data_quick())
    })
    cond.names <- reactive({
        req(filtered_metadata())
        c(
            filtered_metadata() %>% filter(group == unique(filtered_metadata()$group)[unique(filtered_metadata()$group) == selected_conditions()[2]]) %>% pull(intensity_sample_name),
            filtered_metadata() %>% filter(group == unique(filtered_metadata()$group)[unique(filtered_metadata()$group) == selected_conditions()[1]]) %>% pull(intensity_sample_name)
        )
    })

    unique_proteins_reactive <- reactive({
        dq_upf <- data_quick_upf()
        if (input$labeltype == 2) dq_upf <- dq_upf[rowSums(is.na(dq_upf[, cond.names()])) != length(cond.names()), ]
        obtain_unique_proteins(dq_upf, filtered_metadata(), selected_conditions())
    })

    unique_control <- reactive({
        as.data.frame(unique_proteins_reactive()[[1]], check.names = FALSE)
    })
    unique_treatment <- reactive({
        as.data.frame(unique_proteins_reactive()[[2]], check.names = FALSE)
    })

    data_filtered <- reactive({
        df.F <- filter_valids(data_quick_upf(), filtered_metadata(), unique_proteins_reactive(), input$min_prop, at_least_one = FALSE, input$labeltype)
        if (input$uniquefilter) df.F <- unique_peptides_filter(df.F, filtered_metadata(), input$numberuniquepep, input$proportionsamples)
        return(df.F)
    })

    data <- reactive({
        df.F <- data_filtered()
        if (input$normchoice) {
            df.F <- if (input$normoptions == "Median centering") median_centering(df.F, LOG2.names()) else normalization_func(df.F, LOG2.names(), input$normoptions)
        }
        req(input$imputation)
        if (input$imputation == 1) impute_KNN_data(as.data.frame(df.F), LOG2.names(), k = 5) else if (input$imputation == 2) impute_data(as.data.frame(df.F), LOG2.names()) else df.F
    })

    total_dataset <- reactive({
        bind_rows(data(), unique_control()) %>% bind_rows(unique_treatment())
    })

    difexpression <- reactive({
        req(input$pvaladj)
        stat_df <- statistical_analysis(data(), input$displaytest, input$testid, filtered_metadata(), input$LogFC, input$sigcutoff, input$pvaladj, input$statselected, unique_proteins_reactive(), input$proteins, input$PSMaware, input$comptplatform, selected_conditions(), diann_dir = if (input$comptplatform == 3) selected_dir() else NULL)
        stat_df <- stat_df[, !(names(stat_df) %in% "Protein_description")]
        merged <- merge(total_dataset(), stat_df, by = "Protein")
        merged <- if (input$PSMaware) merged %>% select("Protein", "Protein_description", "logFC", "sca.P.Value", "sca.adj.pval", "expression", everything()) else merged %>% select("Protein", "Protein_description", "logFC", "p.value", "adj.P.Val", "expression", everything())
        merged <- merged[order(merged$expression), ]
        row.names(merged) <- merged$Protein
        return(merged)
    })

    # ----------------------------------------
    # SUB-SERVER ROUTERS (Modular Functions)
    # ----------------------------------------
    .server_home(input, output, session)
    .server_data_handling(input, output, session, data, unique_control, unique_treatment, total_dataset, data_quick, LOG2.names, data_filtered, unique_proteins_reactive, filtered_metadata, selected_conditions)
    .server_quality_control(input, output, session, data, filtered_metadata, selected_conditions, LOG2.names, data_filtered)
    .server_differential_analysis(input, output, session, difexpression)
    .server_differential_plots(input, output, session, difexpression, data, filtered_metadata, LOG2.names, selected_conditions)
    .server_functional_analysis(input, output, session, difexpression, data_quick)
    .server_interaction_analysis(input, output, session, difexpression)
    .server_reporting(input, output, session, data_filtered, unique_proteins_reactive, data_quick, filtered_metadata, data, LOG2.names, difexpression, selected_conditions)
}

# ==========================================
# SUB-SERVER DEFINITIONS
# ==========================================

.server_home <- function(input, output, session) {
    output$home_img <- renderImage(
        {
            list(src = system.file("extdata", "traianprot.png", package = "TraianProt"), width = "60%", height = 1000)
        },
        deleteFile = FALSE
    )

    output$download_tutorial <- downloadHandler(
        filename = "TraianProt_Tutorial.pdf",
        content = function(file) {
            file.copy(system.file("extdata", "Tutorial.pdf", package = "TraianProt"), file)
        }
    )
    output$dl_example_matrix <- downloadHandler(
        filename = "proteinGroups.txt",
        content = function(file) {
            file.copy(system.file("extdata", "proteinGroups.txt", package = "TraianProt"), file)
        }
    )
    output$dl_example_metadata <- downloadHandler(
        filename = "metadata_MaxQuant.tsv",
        content = function(file) {
            file.copy(system.file("extdata", "metadata_MaxQuant.tsv", package = "TraianProt"), file)
        }
    )
}

.server_data_handling <- function(input, output, session, data, unique_control, unique_treatment, total_dataset, data_quick, LOG2.names, data_filtered, unique_proteins_reactive, filtered_metadata, selected_conditions) {
    output$file <- DT::renderDataTable({
        req(input$displaytable > 0, input$displayproteins)
        df_disp <- if (input$displayproteins == 1) data() else if (input$displayproteins == 2) unique_control() else unique_treatment()
        DT::datatable(df_disp, options = list(pageLength = 5, lengthMenu = c(5, 10, 15, 20), scrollX = TRUE, autoWidth = TRUE))
    })

    output$downloaddfcomm <- downloadHandler(filename = function() {
        paste(input$filenamedownloadcomm, Sys.Date(), ".txt", sep = "")
    }, content = function(file) {
        write.table(total_dataset(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    })
    output$downloaddfcontroluniques <- downloadHandler(filename = function() {
        paste(input$filenamedownloadcontroluniques, Sys.Date(), ".txt", sep = "")
    }, content = function(file) {
        write.csv(unique_control(), file)
    })
    output$downloaddftreatmentuniques <- downloadHandler(filename = function() {
        paste(input$filenamedownloadtreatmentuniques, Sys.Date(), ".txt", sep = "")
    }, content = function(file) {
        write.csv(unique_treatment(), file)
    })

    output$LOG2.names <- renderTable(
        {
            tryCatch(
                {
                    LOG2.names()
                },
                error = function(e) {
                    data.frame(Message = "Something has gone wrong, check for the settings chosen.")
                }
            )
        },
        striped = TRUE,
        align = "c",
        bordered = TRUE
    )

    proteins_identified <- reactive({
        req(data_quick(), filtered_metadata())
        identify_proteins(data_quick(), filtered_metadata(), input$comptplatform, selected_conditions())
    })

    output$protident <- renderPlot({
        req(p_id <- proteins_identified())
        plot(p_id)
    })

    output$downloadprotquant <- downloadHandler(
        filename = function() "Proteins_quantified.tiff",
        content = function(file) {
            p_id <- proteins_identified()
            grDevices::tiff(file, width = 12, height = 10, units = "in", res = 400)
            plot(p_id)
            grDevices::dev.off()
        }
    )
    power_curve <- reactive({
        traianprot_power_curve(data_quick(), LOG2.names(), input$FoldChangepower, input$replicatespower, input$powerchoice, input$statisticalpower)
    })
    output$powertest <- renderPlot({
        req(power_curve())
        plot(power_curve())
    })

    output$venn <- renderPlot({
        try(grid::grid.draw(venn_diagram(data_filtered(), unique_proteins_reactive(), input$color1, input$color2)), silent = TRUE)
    })
    output$downloadvenn <- downloadHandler(filename = function() {
        paste(input$venn_name, ".", input$vennextension, sep = "")
    }, content = function(file) {
        grDevices::tiff(file, width = 12, height = 10, units = "in", res = 400)
        grid::grid.draw(venn_diagram(data_filtered(), unique_proteins_reactive(), input$color1, input$color2))
        grDevices::dev.off()
    })
}

.server_quality_control <- function(input, output, session, data, filtered_metadata, selected_conditions, LOG2.names, data_filtered) {
    boxplot_distribution <- reactive({
        boxplot_function(data(), filtered_metadata(), selected_conditions())
    })
    pre_imp_plot <- reactive({
        preimputation_state(data_filtered(), filtered_metadata()$log2_col)
    })
    correlation_plot <- reactive({
        corrplot_function(data()[filtered_metadata()$log2_col], filtered_metadata(), input$dispmethod)
    })

    output$boxplot <- renderPlot({
        if (input$displaydistplots == 1) {
            try(plot(boxplot_distribution()), silent = TRUE)
        } else if (input$displaydistplots == 2) try(plotCV2(data()[, LOG2.names()], trend = TRUE, main = "Dispersion check", cex = 0.2, pch = 16, xlab = "Average log-intensity", ylab = expression("Relative standard deviation")), silent = TRUE)
    })

    output$preimputationplot <- renderPlot({
        if (input$displayimpplots == 1) {
            try(plot(pre_imp_plot()), silent = TRUE)
        } else if (input$displayimpplots == 2) try(plot(postimputation_state(data_filtered(), input$imputation, filtered_metadata()$log2_col)), silent = TRUE)
    })

    output$pcaplot <- renderPlot({
        if (input$dimreduction == 1) {
            try(plot(pca(data(), filtered_metadata(), selected_conditions(), pc_x = input$pc_x_axis, pc_y = input$pc_y_axis)), silent = TRUE)
        } else if (input$dimreduction == 2) try(plot(tsne(data(), filtered_metadata(), perplexity_num = input$perplexitytsne, selected_conditions())), silent = TRUE)
    })

    output$histogram <- renderPlot({
        if (input$displaynormplots == 1) {
            histogram(data(), input$histsample, input$histcol, input$histtitle)
        } else if (input$displaynormplots == 2) qqplot_function(data(), input$scatsample1, input$scatsample2, input$qqcolor)
    })

    output$sscatplot <- renderPlot({
        if (input$displaycorrplots == 1) {
            scatterplot_function(data(), input$scatsample1, input$scatsample2)
        } else if (input$displaycorrplots == 2) try(plot(correlation_plot()), silent = TRUE)
    })

    output$downloadqmplot <- downloadHandler(
        filename = function() paste0(input$qualitymetricsfilename, ".tiff"),
        content = function(file) {
            choice <- as.character(input$qualitymetricsplotchoice)

            r <- if (choice == "2") 300 else 400
            w_px <- if (choice == "2") 1500L else 4800L
            h_px <- if (choice == "2") 1200L else 4000L

            plot_obj <- if (choice != "2") {
                switch(choice,
                    "1" = boxplot_function(data(), filtered_metadata(), selected_conditions()),
                    "3" = preimputation_state(data_filtered(), filtered_metadata()$log2_col),
                    "4" = postimputation_state(data_filtered(), input$imputation, filtered_metadata()$log2_col),
                    "5" = histogram(data(), input$histsample, input$histcol, input$histtitle),
                    "6" = qqplot_function(data(), input$scatsample1, input$scatsample2, input$qqcolor),
                    "7" = scatterplot_function(data(), input$scatsample1, input$scatsample2),
                    "8" = corrplot_function(data()[filtered_metadata()$log2_col], filtered_metadata(), input$dispmethod),
                    "9" = pca(data(), filtered_metadata(), selected_conditions(), pc_x = input$pc_x_axis, pc_y = input$pc_y_axis),
                    "10" = tsne(data(), filtered_metadata(), perplexity_num = input$perplexitytsne, selected_conditions())
                )
            }
            if (inherits(plot_obj, "ggplot")) {

                ggplot2::ggsave(file,
                    plot = plot_obj, device = "tiff", width = w_px, height = h_px,
                    units = "px", dpi = r, limitsize = FALSE, compression = "lzw"
                )
            } else {
                grDevices::tiff(file, width = w_px, height = h_px, units = "px", res = r, compression = "lzw")
                if (choice == "2") {
                    plotCV2(data()[, LOG2.names()],
                        trend = TRUE, main = "Dispersion check", cex = 0.2,
                        pch = 16, xlab = "Average log-intensity", ylab = expression("Relative standard deviation")
                    )
                } else {
                    plot(plot_obj)
                }
                grDevices::dev.off()
            }
        }
    )
}

.server_differential_analysis <- function(input, output, session, difexpression) {
    output$significantBox_dm <- renderInfoBox({
        num_total <- nrow(difexpression())
        num_signif <- nrow(subset(difexpression(), difexpression()$expression != "Unchanged"))
        frac <- num_signif / num_total
        infoBox("Significant proteins", paste0(num_signif, " out of ", num_total), paste0(signif(frac * 100, digits = 3), "% of proteins differentially expressed"), icon = icon("stats", lib = "glyphicon"), color = "blue", width = 10)
    })

    output$limma <- DT::renderDataTable({
        DT::datatable(difexpression(), options = list(pageLength = 15, lengthMenu = c(5, 10, 15, 20), scrollX = TRUE, autoWidth = TRUE))
    })
    output$downloaddif <- downloadHandler(filename = function() {
        paste(input$namedownload, Sys.Date(), ".csv", sep = "")
    }, content = function(file) {
        write.csv(difexpression(), file)
    })
}

.server_differential_plots <- function(input, output, session, difexpression, data, filtered_metadata, LOG2.names, selected_conditions) {
    output$volcanoplot <- renderPlotly({
        volcano_plot(difexpression(), input$volcanotitle, input$labelprots, input$statselected, input$PSMaware)
    })
    output$heatmapplot <- renderPlot({
        tryCatch(
            {
                my_heatmap(data(), filtered_metadata()$log2_col, input$heatmaptitle)
            },
            error = function(e) {
                plot.new()
                text(0.5, 0.5, "Error: check missing values density.", cex = 1.2)
            }
        )
    })
    output$difheatmapplot <- renderPlot({
        tryCatch(
            {
                (my_heatmap_differential(difexpression(), data(), filtered_metadata()$log2_col, input$heatmaptitle))
            },
            error = function(e) {
                plot.new()
                text(0.5, 0.5, "Amount of NAs too large.", cex = 1.2)
            }
        )
    })
    output$difboxplotplot <- renderPlot({
        Diferential_boxplot(data(), filtered_metadata(), protein = input$proteinname, LOG2.names(), selected_conditions())
    })

    output$downloadvolcano <- downloadHandler(
        filename = function() {
            paste0(input$name_download_volcano, ".", input$difextension)
        },
        content = function(file) {
            p_volc <- volcano_plot_tiff(difexpression(), input$volcanotitle, input$labelprots, input$statselected, input$PSMaware)
            do.call(getExportedValue("grDevices", input$difextension), c(list(file, width = 12, height = 10), if (input$difextension != "pdf") list(units = "in", res = 400)))
            plot(p_volc)
            grDevices::dev.off()
        }
    )

    output$downloadheatmap <- downloadHandler(
        filename = function() {
            paste0(input$name_download_heatmap, ".", input$difextension)
        },
        content = function(file) {
            p_heat <- my_heatmap(data(), filtered_metadata()$log2_col, input$heatmaptitle)

            do.call(getExportedValue("grDevices", input$difextension), c(list(file, width = 12, height = 10), if (input$difextension != "pdf") list(units = "in", res = 400)))
            methods::show(p_heat)
            grDevices::dev.off()
        }
    )

    output$downloaddifheatmap <- downloadHandler(
        filename = function() {
            paste0(input$name_download_heatmap, ".", input$difextension)
        },
        content = function(file) {
            p_heat <- my_heatmap_differential(difexpression(), data(), filtered_metadata()$log2_col, input$heatmaptitle)
            do.call(getExportedValue("grDevices", input$difextension), c(list(file, width = 12, height = 10), if (input$difextension != "pdf") list(units = "in", res = 400)))
            methods::show(p_heat)
            grDevices::dev.off()
        }
    )
}

.server_functional_analysis <- function(input, output, session, difexpression, data_quick) {
    funcanalysis <- reactive({
        req(input$organism, input$proteinid)
        Goterms_finder(difexpression(), data_quick(), target = input$proteinid, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organismo = input$organism, custombg = input$backgroundset, input$comptplatform, user_threshold = input$userthreshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC"))
    })

    output$functionalanalysis <- DT::renderDataTable({
        req(funcanalysis())
        DT::datatable(funcanalysis()[[2]]@result, options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE))
    })
    output$downloadfunc <- downloadHandler(filename = function() {
        paste0(input$name_download, Sys.Date(), ".xlsx")
    }, content = function(file) {
        writexl::write_xlsx(funcanalysis()[[2]]@result, file)
    })

    dotplot_react <- reactive({
        req(funcanalysis())
        dotplot_func(funcanalysis(), x = "GeneRatio", title = input$dotplottitle, split = "Conditions", font.size = input$fontsizedotplot, showCategory = input$categorydotplot, color = "adj.P.Val")
    })

    output$dotplotout <- renderPlot({
        req(p_dot <- dotplot_react())
        p_dot
    })

    output$downloaddotplot <- downloadHandler(
        filename = function() {
            paste0(input$name_download_func, ".", input$funcextension)
        },
        content = function(file) {
            quality_input <- input$funcquality
            dpi_val <- if (!is.null(quality_input) && grepl("^\\d+$", quality_input)) as.numeric(quality_input) else 300

            ggplot2::ggsave(file,
                plot = dotplot_react(), device = input$funcextension,
                width = 12, height = 10, units = "in", dpi = dpi_val
            )
        }
    )
    barplot_react <- reactive({
        req(funcanalysis())
        barplot_func(funcanalysis(), input$categorybarplot, conditions = input$barplottitle, font.size = input$fontsizebarplot)
    })

    output$barplotout <- renderPlot({
        req(p_bar <- barplot_react())
        p_bar
    })

    output$downloadbarplot <- downloadHandler(
        filename = function() {
            paste0(input$name_download_func, ".", input$funcextension)
        },
        content = function(file) {
            quality_input <- input$funcquality
            dpi_val <- if (!is.null(quality_input) && grepl("^\\d+$", quality_input)) as.numeric(quality_input) else 300

            ggplot2::ggsave(file,
                plot = barplot_react(), device = input$funcextension,
                width = 12, height = 10, units = "in", dpi = dpi_val
            )
        }
    )

    output$manhattanout <- renderPlotly({
        req(input$manhattan > 0)
        gostplot_func(funcanalysis())
    })
}

.server_interaction_analysis <- function(input, output, session, difexpression) {
    interactions <- reactive({
        interactions_up(difexpression(), input$taxonid, input$scthreshold)
    })
    output$upreg <- renderPlot({
        req(input$upnet > 0)
        interactions_up(difexpression(), input$taxonid, input$scthreshold)
    })
    output$downreg <- renderPlot({
        req(input$downnet > 0)
        interactions_down(difexpression(), input$taxonid, input$scthreshold)
    })

    graph <- reactive({
        igraph_analysis(interactions(), input$taxonid, input$scthreshold)
    })
    output$upgraph <- DT::renderDataTable({
        datatable(graph()[[1]], options = list(pageLength = 1, scrollX = TRUE))
    })
    output$downgraph <- DT::renderDataTable({
        datatable(graph()[[2]], options = list(pageLength = 1, scrollX = TRUE))
    })
}

.server_reporting <- function(input, output, session, data_filtered, unique_proteins_reactive, data_quick, filtered_metadata, data, LOG2.names, difexpression, selected_conditions) {
    output$downloadreport <- downloadHandler(
        filename = function() {
            paste0("Downstreaming_Analysis_report-", Sys.Date(), ".html")
        },
        content = function(file) {
            tryCatch(
                {
                    tempReport <- file.path(tempdir(), "Downstreaming_Analysis_report.Rmd")
                    file.copy(system.file("extdata", "Downstreaming_Analysis_report.Rmd", package = "TraianProt"), tempReport, overwrite = TRUE)
                    fa <- Goterms_finder(difexpression(), data_quick(), target = input$proteinid, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organismo = input$organism, custombg = input$backgroundset, input$comptplatform, user_threshold = input$userthreshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC"))

                    params <- list(
                        venn = venn_diagram(data_filtered(), unique_proteins_reactive(), input$color1, input$color2),
                        proteins = identify_proteins(data_quick(), filtered_metadata(), input$comptplatform, selected_conditions()),
                        boxplot = boxplot_function(data(), filtered_metadata(), selected_conditions()),
                        pca = pca(data(), filtered_metadata(), selected_conditions(), pc_x = input$pc_x_axis, pc_y = input$pc_y_axis), # <- Arreglado también aquí
                        imputation = preimputation_state(data_filtered(), filtered_metadata()$log2_col),
                        corrplot = corrplot_function(data()[filtered_metadata()$log2_col], filtered_metadata(), input$dispmethod),
                        volcano = volcano_plot(difexpression(), input$volcanotitle, input$labelprots, input$statselected, input$PSMaware),
                        heatmap = my_heatmap_differential(difexpression(), data(), filtered_metadata()$log2_col, input$heatmaptitle),
                        dotplot = dotplot_func(fa, x = "GeneRatio", title = input$dotplottitle, split = "Conditions", font.size = input$fontsizedotplot, showCategory = input$categorydotplot, color = "adj.P.Val"),
                        barplot = barplot_func(fa, input$categorybarplot, conditions = input$barplottitle, font.size = input$fontsizebarplot),
                        control_replicates = input$repcond1,
                        treatment_replicates = input$repcond2,
                        proportion_filtering_int = input$min_prop,
                        proportion_filtering_unique_peptides = input$proportionsamples,
                        unique_peptides = input$numberuniquepep,
                        type_of_test = input$testid,
                        log2FC = input$LogFC,
                        pvalue = input$sigcutoff,
                        organism = input$organism,
                        p_value_for_enrichment = input$userthreshold,
                        taxon_id = input$taxonid,
                        score_threshold = input$scthreshold
                    )

                    rmarkdown::render(tempReport, output_format = "html_document", output_file = file, params = params, envir = new.env(parent = globalenv()))
                },
                error = function(e) {
                    message("Something went wrong: ", e$message)
                }
            )
        }
    )
}
