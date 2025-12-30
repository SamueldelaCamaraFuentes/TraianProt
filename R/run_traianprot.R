#' Ejecuta el pipeline completo de TraianProt
#' @export
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
    logfcup = 1,
    logfcdown = -1,
    sig = 0.05,
    conditions = NULL,  # Si es NULL, se genera automático
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
  
  # --- GESTIÓN DE TÍTULO AUTOMÁTICO ---
  if (is.null(conditions)) {
    conditions <- paste(selected_conditions, collapse = " vs ")
  }
  
  # 1. CARGA DE DATOS
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
  
  # 2. Selección de Condiciones (Lógica de filtrado de grupos)
  
  # Listar condiciones disponibles
  available_conditions <- unique(metadata$group)
  message("Condiciones encontradas en metadata: ", paste(available_conditions, collapse = ", "))
  
  # Verificar que las condiciones existen
  if(!all(selected_conditions %in% available_conditions)) {
    stop("Alguna de las condiciones seleccionadas no existe en la columna 'Group' del metadata.")
  }
  
  # Filtrar el metadata
  filtered_metadata <- metadata[metadata$group %in% selected_conditions, ]
  message("Metadata filtrado. Muestras seleccionadas: ", nrow(filtered_metadata))
  
  # Carga de RAW Data según plataforma
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
  
  # 4.2) Quantification filtering (Corregido: at_least_one = FALSE)
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
  
  # Histogram (Solo imprime en pantalla/RStudio, no guarda tiff si no se envuelve)
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
  # CORREGIDO: Eliminada linea colgada y arreglado nombre
  filename_tsne <- file.path(output_directory, "TSNE.tiff")
  tiff(filename = filename_tsne, units="in", width=12, height=10, res=400)
  tsne(df.FNI, filtered_metadata, perplexity_num = perplexity_num, selected_conditions)
  dev.off()
  
  
  # 6. Differential analysis
  # ------------------------------------------------------------------------------
  pvaladj <- "fdr"
  
  limma_1 <- statistical_analysis(df.FNI, test, paired, filtered_metadata, logfcup, logfcdown, sig, pvaladj, statval,
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
  
  # CORREGIDO: Guardado dinámico del CSV con resultados estadísticos
  csv_name <- paste0("Differential_Analysis_", selected_conditions[1], "_vs_", selected_conditions[2], ".csv")
  write.csv(limma, file = file.path(output_directory, csv_name), row.names = FALSE)
  
  
  # Representaciones 
  
  # Volcano
  volcano <- volcano_plot_tiff(limma, conditions, label, statval, psms)
  filename <- file.path(output_directory, "Volcano.tiff")
  
  tiff(filename = filename, units="in", width=12, height=10, res=400) 
  print(volcano)
  dev.off()
  
  # Heatmap
  # CORREGIDO: Título dinámico
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
  
  # Representaciones
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
  
  # Si deseas guardar el grafo, deberías añadir código aquí para plotearlo o guardarlo.
  graph_analysis <- igraph_analysis(interactions, taxonid, score)
  
  message("Pipeline finalizado con éxito.")
}





