#' Run the Complete TraianProt Analysis Pipeline
#'
#' @description
#' This function executes the full TraianProt workflow for proteomics data analysis.
#' It performs data loading, filtering, normalization, missing value imputation,
#' differential expression analysis, functional enrichment analysis (GO/KEGG),
#' and protein-protein interaction (PPI) network analysis. Results are saved
#' as Excel files and high-resolution plots (TIFF) in the specified output directory.
#'
#' @details
#' The pipeline consists of the following steps:
#' \enumerate{
#'   \item \strong{Data Loading & Preprocessing:} Imports raw data and metadata.
#'   \item \strong{Filtering:} Removes contaminants, reverse hits, and proteins with low valid values.
#'   \item \strong{Normalization:} Applies normalization methods (e.g., "mean", "median") to reduce technical bias.
#'   \item \strong{Imputation:} Imputes missing values using methods like KNN or MinProb.
#'   \item \strong{Differential Analysis:} Uses \code{limma} or \code{DEqMS} to identify significantly changed proteins.
#'   \item \strong{Functional Enrichment:} Performs GO/KEGG enrichment analysis using \code{gprofiler2}.
#'   \item \strong{PPI Analysis:} Retrieves interaction networks from STRINGdb.
#' }
#'
#' @param file_path_raw Character string. Path to the raw proteomics data file (e.g., "proteinGroups.txt").
#' @param file_path_metadata Character string. Path to the metadata file (CSV/XLSX) containing "intensity_sample_name", "group", and "log2_col".
#' @param output_directory Character string. Path to the folder where results and plots will be saved.
#' @param selected_conditions Character vector. The two experimental conditions to compare (e.g., \code{c("Control", "Treatment")}).
#' @param file_path_report Character string. Path to the report file (only required for specific platforms like Proteome Discoverer). Default is \code{NULL}.
#' @param platform_id Numeric. Identifier for the software used to generate raw data:
#' \itemize{
#'   \item 1: MaxQuant
#'   \item 2: MSFragger
#'   \item 3: Proteome Discoverer
#'   \item 4: DIA-NN
#'   \item 5: Spectronaut/ProteoScape
#' }
#' @param organism_id Numeric. Identifier for the organism (used for filtering contaminants): 1 (Candida albicans), 2 (Other), 3 (Rat).
#' @param labeltype Numeric. Type of labeling: 1 (Label-Free), 2 (TMT/Label-based).
#' @param norm_method Character string. Normalization method to apply. Options: "mean", "median", "quantile", or "none". Default is "mean".
#' @param imputation_method Character string. Missing value imputation method:
#' \itemize{
#'   \item "1": MinProb (probabilistic minimum)
#'   \item "2": QRILC (Quantile Regression Imputation)
#'   \item "3": KNN (K-Nearest Neighbors) - Default
#' }
#' @param min_prop_filtering Numeric. Minimum proportion of valid values required per group (0 to 1). Default is 0.5.
#' @param min_unique_peptides Numeric. Minimum number of unique peptides required to keep a protein. Default is 1.
#' @param test Numeric. Statistical test to perform: 1 (T-test), 2 (Limma) - Default.
#' @param psms Logical. Whether to use PSM count correction (DEqMS) if available. Default is \code{TRUE}.
#' @param paired Logical. Perform paired analysis? Default is \code{FALSE}.
#' @param statval Numeric. Statistic used for thresholding: 1 (P-value), 2 (Adjusted P-value). Default is 1.
#' @param way Numeric. Direction of analysis: 1 (One-way), 2 (Two-way). Default is 2.
#' @param logfc Numeric. Log2 Fold Change threshold for upregulation and downregulation (e.g., 0.58 for 1.5-fold). Default is 1.
#' @param sig Numeric. Significance threshold (alpha). Default is 0.05.
#' @param conditions Character string (Optional). Custom label for the comparison (e.g., "Treated vs Control"). If \code{NULL}, it is auto-generated.
#' @param id_target Character string. Gene ID type for enrichment (e.g., "ENSG", "ENTREZGENE"). Default is "ENSG".
#' @param organismo Character string. Organism name for \code{gprofiler2} (e.g., "hsapiens", "mmusculus").
#' @param threshold Numeric. Significance threshold for enrichment analysis. Default is 0.05.
#' @param custombg Logical. Whether to use the identified proteins as a custom background for enrichment. Default is \code{FALSE}.
#' @param label Numeric. Labeling style for plots. Default is 3.
#' @param font_size Numeric. Font size for generated plots. Default is 12.
#' @param terms_number Numeric. Number of top enrichment terms to display in plots. Default is 12.
#' @param taxonid Numeric. NCBI Taxonomy ID for STRINGdb (e.g., 9606 for Human). Default is 9606.
#' @param score Numeric. Confidence score threshold for STRING interactions (0-1000). Default is 400.
#'
#' @return No return value. The function saves the following files in \code{output_directory}:
#' \itemize{
#'   \item \strong{Plots:} PCA, Volcano Plot, Heatmaps, Boxplots, Correlation plots (TIFF format).
#'   \item \strong{Tables:} Differential expression results, Enrichment results (Excel format).
#' }
#'
#' @importFrom grDevices tiff dev.off
#' @importFrom graphics barplot
#' @importFrom utils write.csv
#' @importFrom stats complete.cases median
#' @export
#'
#' @examples
#' \donttest{
#' # Example usage assuming files exist in the working directory
#' run_traianprot(
#'   file_path_raw = "data/proteinGroups.txt",
#'   file_path_metadata = "data/metadata.xlsx",
#'   output_directory = "results/",
#'   selected_conditions = c("Control", "Treatment"),
#'   platform_id = 1, # MaxQuant
#'   organism_id = 1, # Human
#'   imputation_method = "3", # KNN
#'   organismo = "hsapiens"
#' )
#' }
run_traianprot <- function(
    file_path_raw,
    file_path_metadata,
    output_directory,
    selected_conditions, # Example: c("Treatment", "Control")

    file_path_report = NULL, # Only required for platform 3
    platform_id = 1,
    organism_id = 1,
    labeltype = 1,
    norm_method = "mean",
    imputation_method = "3",
    min_prop_filtering = 0.5,
    min_unique_peptides = 1,
    test = 2,
    psms = TRUE,
    paired = FALSE,
    statval = 1,
    way = 2,
    logfc = 1,
    sig = 0.05,
    conditions = NULL,
    id_target = "ENSG",
    organismo = NULL,
    threshold = 0.05,
    custombg = FALSE,
    label = 3,
    font_size = 12,
    terms_number = 12,
    taxonid = 9606, # Default: Human
    score = 400
) {


  if (is.null(conditions)) {
    conditions <- paste(selected_conditions, collapse = " vs ")
  }

  # 1. Data Loading
  # ------------------------------------------------------------------------------

  # Metadata
  metadata <- read.delim(file_path_metadata, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

  if (platform_id == 1 | platform_id == 2 | platform_id == 5){
    metadata <- metadata %>%
      mutate(
        raw_name = intensity_sample_name,
        intensity_sample_name = raw_name,
        log2_col = sub("Intensity", "LOG2", raw_name)
      ) %>%
      mutate(
        raw_name = intensity_sample_name,
        unique_peptides_col = sub("Intensity", "Unique.peptides", raw_name)
      )%>%
      select(intensity_sample_name, group, sample_name, log2_col, unique_peptides_col)

  } else if (platform_id == 3) {
    metadata <- metadata %>%
      mutate(
        raw_name = basename(intensity_sample_name),
        intensity_sample_name = raw_name,
        log2_col = sub("\\.(d|raw)$", ".LOG2", raw_name, ignore.case = TRUE),
        unique_peptides_col = paste0("Unique peptides ", sub("\\.(d|raw)$", "", raw_name, ignore.case = TRUE))
      ) %>%
      select(intensity_sample_name, group, sample_name, log2_col, unique_peptides_col)

  } else if (platform_id == 4){
    metadata <- metadata %>%
      mutate(
        raw_name = intensity_sample_name,
        intensity_sample_name = raw_name,
        log2_col = sub("Abundance:", "LOG2", raw_name)
      ) %>%
      mutate(
        raw_name = intensity_sample_name,
        unique_peptides_col = sub("Abundance:", "Unique.peptides", raw_name)
      )%>%
      select(intensity_sample_name, group, sample_name, log2_col, unique_peptides_col)
  }

  # 2. Conditions selection (Filtering group logic)

  available_conditions <- unique(metadata$group)
  message("Condiciones encontradas en metadata: ", paste(available_conditions, collapse = ", "))

  # Verification
  if(!all(selected_conditions %in% available_conditions)) {
    stop("Alguna de las condiciones seleccionadas no existe en la columna 'Group' del metadata.")
  }

  filtered_metadata <- metadata[metadata$group %in% selected_conditions, ]
  message("Metadata filtrado. Muestras seleccionadas: ", nrow(filtered_metadata))

  if (platform_id == 1 | platform_id == 2 | platform_id == 5){
    raw <- read.delim(file_path_raw, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
    df <- quick_filtering(raw, platform_id, organism_id, filtered_metadata, selected_conditions)

  } else if (platform_id == 3){
    raw <- read.delim(file_path_raw, sep = "\t", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
    df <- quick_filtering(raw, platform_id, organism_id, filtered_metadata, selected_conditions, file_path_report)
    df <- as.data.frame(df)

  } else if (platform_id == 4){
    raw <- as.data.frame(readxl::read_xlsx(file_path_raw))
    df <- quick_filtering(raw, platform_id, organism_id, filtered_metadata, selected_conditions)
  }

  # 3. LOG 2 INTENSITY
  LOG2.names <- obtain_LOG.names(df)

  # Unique proteins
  unique_proteins <- obtain_unique_proteins(df, filtered_metadata, selected_conditions)

  # 4. Preprocessing
  # ------------------------------------------------------------------------------

  # 4.1) UNIQUE PEPTIDES FILTERING
  df <- unique_peptides_filter(df, filtered_metadata, min_unique_peptides, min_prop_filtering)

  # 4.2) Quantification filtering
  df.F <- filter_valids(df, filtered_metadata, unique_proteins, min_prop_filtering, at_least_one = FALSE, labeltype)

  # 4.3) Normalization
  df.F <- median_centering(df.F, LOG2.names)
  df.F <- normalization_func(df.F, LOG2.names, norm_method)

  # 4.4) Imputation
  df.FNI <- impute_KNN_data(as.data.frame(df.F), LOG2.names, k = 5)
  df.FNI <- impute_data(as.data.frame(df.F), LOG2.names)

  total <- bind_rows(df.FNI, as.data.frame(unique_proteins[1], check.names = FALSE))
  total_dataset <- bind_rows(total, as.data.frame(unique_proteins[2], check.names = FALSE))


  ###############
  # Venn
  ##############
  filename <- file.path(output_directory, "Venn Diagram.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  grid.draw(venn_diagram(df.F, unique_proteins, label1 = "Control", label2 = "KO", color1 = "blue", color2 = "maroon"))
  dev.off()

  #####################
  # Proteins identified
  #####################
  filename <- file.path(output_directory, "Proteins_identified.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  identify_proteins(df, filtered_metadata, platform_id, selected_conditions)
  dev.off()

  # 5. Quality metrics
  # ------------------------------------------------------------------------------

  filename <- file.path(output_directory, "CV2plot.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  plotCV2(df.FNI[,LOG2.names], trend = TRUE, main = "Dispersion check", cex = 0.2, pch = 16, xlab="Average log-intensity", ylab=expression("Relative standard deviation"))
  dev.off()

  # Boxplot
  filename <- file.path(output_directory, "Boxplot.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  boxplot <- boxplot_function(df.FNI, filtered_metadata, selected_conditions)
  print(boxplot)
  dev.off()

  # Check NAs
  filename <- file.path(output_directory, "Preimputation_state.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  preimputation_state(df.F, filtered_metadata$log2_col)
  dev.off()

  filename <- file.path(output_directory, "Postimputation_state.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  preimputation_state(df.FNI, filtered_metadata$log2_col)
  dev.off()

  # Histogram
  colname <- filtered_metadata$log2_col[1]
  color <- "maroon"
  title <- "Histogram of intensities"
  histogram(df.FNI, colname, color, title)

  # Scatter plot
  sample1 <- filtered_metadata$log2_col[1]
  sample2 <- filtered_metadata$log2_col[2]

  filename1 <- file.path(output_directory, "Scatterplot_definitivo.tiff")
  tiff(filename = filename1, units="in", width=12, height=10, res=400)
  scatterplot_function(df.FNI, sample1, sample2)
  dev.off()

  # QQ plot
  filename2 <- file.path(output_directory, "QQplot.tiff")
  tiff(filename = filename2, units="in", width=12, height=10, res=400)
  qq_plot <- qqplot_function(df.FNI, sample1, sample2, color = "blue")
  print(qq_plot)
  dev.off()

  # Correlation plot
  filename3 <- file.path(output_directory, "Correlation_plot.tiff")
  tiff(filename = filename3, units="in", width=12, height=10, res=400)
  corrplot_function(df.FNI[filtered_metadata$log2_col], display = "shade")
  dev.off()

  # PCA
  filename3 <- file.path(output_directory, "PCA.tiff")
  tiff(filename = filename3, units="in", width=12, height=10, res=400)
  pca(df.FNI,filtered_metadata, selected_conditions, pc_x = 1, pc_y = 2)
  dev.off()

  # T-SNE
  perplexity_num <- 2

  filename_tsne <- file.path(output_directory, "TSNE.tiff")
  tiff(filename = filename_tsne, units="in", width=12, height=10, res=400)
  tsne(df.FNI, filtered_metadata, perplexity_num = perplexity_num, selected_conditions)
  dev.off()


  # 6. Differential analysis
  # ------------------------------------------------------------------------------
  pvaladj <- "fdr"

  limma_1 <- statistical_analysis(df.FNI, test, paired, filtered_metadata, logfc, sig, pvaladj, statval,
                                  unique_proteins, way, psms, platform_id, selected_conditions, diann_dir = if (platform_id == 3) file_path_report else NULL)

  limma <- limma_1[-6]
  limma <- merge(total_dataset, limma, by = "Protein", check.names = FALSE)

  if (psms == TRUE){
    limma <- limma %>%
      dplyr::select("Protein", "Protein_description", "logFC", "sca.P.Value", "sca.adj.pval", "expression", everything())
  } else if (psms == FALSE){
    limma <- limma %>%
      select("Protein", "Protein_description", "logFC", "p.value", "adj.P.Val", "expression", everything())
  }

  limma <- limma[order(limma$expression),]
  row.names(limma) <- limma$Protein

  csv_name <- paste0("Differential_Analysis_", selected_conditions[1], "_vs_", selected_conditions[2], ".csv")
  write.csv(limma, file = file.path(output_directory, csv_name), row.names = FALSE)


  # Differential analysis plots

  # Volcano
  volcano <- volcano_plot_tiff(limma, conditions, label, statval, psms)
  filename <- file.path(output_directory, "Volcano.tiff")

  tiff(filename = filename, units="in", width=12, height=10, res=400)
  print(volcano)
  dev.off()

  # Heatmap
  heatmap_title <- conditions
  heatmap <- my_heatmap(df.FNI, filtered_metadata$log2_col, heatmap_title)

  filename <- file.path(output_directory, "Heatmap.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  my_heatmap(df.FNI, filtered_metadata$log2_col, heatmap_title)
  dev.off()

  diff_heatmap <- my_heatmap_differential(limma, df.FNI, filtered_metadata$log2_col, heatmap_title)
  filename <- file.path(output_directory, "Diferential_Heatmap.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  print(diff_heatmap)
  dev.off()

  # 7. Functional analysis
  # ------------------------------------------------------------------------------

  Go_terms <- Goterms_finder(limma, df, target = id_target, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organismo = organismo, custombg = custombg, platform_id,  user_threshold = threshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC"))

  out_name <- paste0("Enrichment_", selected_conditions[1], "_vs_", selected_conditions[2], ".xlsx")
  writexl::write_xlsx(Go_terms[[2]]@result, path = file.path(output_directory, out_name))

  # Plots
  filename <- file.path(output_directory, "Dotplot.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  dotplot_func(Go_terms, x = "GeneRatio", title = conditions, split = "Conditions", font.size = font_size, showCategory = terms_number, color = "adj.P.Val")
  dev.off()

  filename <- file.path(output_directory, "Barplot.tiff")
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  barplot_func(Go_terms, terms_number, conditions = conditions, font.size = font_size)
  dev.off()

  # 8. Interaction analysis
  # ------------------------------------------------------------------------------

  interactions <- interactions_up(limma, taxonid, score)
  interactions_down <- interactions_down(limma, taxonid, score)

  graph_analysis <- igraph_analysis(interactions, taxonid, score)

  message("Pipeline finalizado con exito.")
}





