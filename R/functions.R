
#################################### helper functions ####################################


##################################################################################
#Quick filtering

#' @title Quick Data Filtering
#' @description Filters raw proteomics data (e.g., from MaxQuant) to remove
#' contaminants, reverse hits, and site-only identifications. It also parses
#' protein identifiers from Fasta headers.
#'
#' @param raw A data frame of raw proteomics data (e.g., from MaxQuant).
#' @param platform An indicator (e.g., 1) for the data platform (e.g., MaxQuant).
#' @param organism An indicator (e.g., 1 for Human, 2 for Mouse, etc.).
#' @param metadata A data frame mapping intensity column names to sample names.
#' @param selected_conditions (Not currently used in function) A vector of
#' selected conditions.
#' @param directory_path (Not currently used in function) A path to a directory.
#'
#' @return A data frame filtered to remove contaminants/reverse hits, with
#' new "Protein" and "Protein_description" columns.
#' @export
#' @importFrom dplyr filter
#' @importFrom utils read.delim
#'
#' @examples
#' # Create synthetic raw data
#' raw_data <- data.frame(
#'   Potential.contaminant = c("", "+", ""),
#'   Reverse = c("", "", ""),
#'   Only.identified.by.site = c("", "", "+"),
#'   Fasta.headers = c(">sp|P12345|GENE1_HUMAN",
#'                     ">sp|P23456|GENE2_HUMAN",
#'                     ">sp|P34567|GENE3_HUMAN"),
#'   Intensity.Sample1 = c(1,2,3),
#'   stringsAsFactors = FALSE
#' )
#' metadata <- data.frame(intensity_sample_name = "Intensity.Sample1")
#'
#' # Run the filtering
#' filtered_df <- quick_filtering(raw_data, 1, 1, metadata, NULL, NULL)
#' print(filtered_df)
#'

quick_filtering <- function(raw, platform, organism, metadata, selected_conditions, directory_path = NULL){
  intensity_names <- metadata$intensity_sample_name
  if (platform == 1){
    df <- raw %>%
      dplyr::filter(Potential.contaminant != "+") %>%
      dplyr::filter(Reverse != "+") %>%
      dplyr::filter(Only.identified.by.site != "+")

    if (organism == 1){
      #Obtenemos los identificadores
      regex <- regexpr(".*(?=.CGDID)", df$Fasta.headers, perl = TRUE)
      df$Protein <- regmatches(df$Fasta.headers, regex)

      df$Protein <- sub(" ", ";", df$Protein)
      regex2 <- regexpr("(?<=;).*(?>[A-Z0-9])", df$Protein, perl = TRUE)
      df$Protein <- regmatches(df$Protein, regex2)

      #Obtenemos la descripci?n proteica
      regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Fasta.headers, perl = TRUE)
      df$Protein_description <- regmatches(df$Fasta.headers, regex3)
    } else if (organism == 2){

      pattern <- regexpr("(?<=\\|)[A-Z0-9-]+(?=\\|)", df$Fasta.headers, perl = TRUE)
      df$Protein <- regmatches(df$Fasta.headers, pattern)

      match <- regexpr("(?<=\\s)(.*?)(?=\\sOS)", df$Fasta.headers, perl = TRUE)
      df$Protein_description <- regmatches( df$Fasta.headers, match)
    }

    #Hacemos la log-transformaci?n

    df[intensity_names] <- sapply(df[intensity_names], as.numeric)
    LOG_names <- sub("Intensity", "LOG2", intensity_names)
    df[LOG_names] <- log2(df[intensity_names])

    return(df)

  } else if (platform == 2){
    df <- raw

    if (organism == 1){
      regex <- regexpr(".*(?=.CGDID)", df$Protein.ID, perl = TRUE)
      df$Protein <- regmatches(df$Protein.ID, regex)

      df$Protein <- sub(" ", ";", df$Protein)
      regex2 <- regexpr("(?<=;).*(?>[A-Z0-9])", df$Protein, perl = TRUE)
      df$Protein <- regmatches(df$Protein, regex2)

      #Obtenemos la descripci?n proteica
      regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Protein.ID, perl = TRUE)
      df$Protein_description <- regmatches(df$Protein.ID, regex3)

      colnames(df)[1] <- "Protein"

    } else if (organism == 2){

      colnames(df)[1] <- "Protein.ID"
      colnames(df)[2] <- "Protein"
      colnames(df)[3] <- "Protein_description"
    }

    df[intensity_names] <- sapply(df[intensity_names], as.numeric)
    LOG2.names <- metadata$log2_col
    df[LOG2.names] <- log2(df[intensity_names])

    return(df)
  } else if (platform == 3){
    df <- raw
    colnames(df)[1] <- "Protein"
    colnames(df)[3] <- "Protein_description"

    intensity_columns <- grepl("\\.(d|raw)$", colnames(df))
    intensity_columns <- colnames(df[intensity_columns])

    #Le damos un lavado
    new_names <- c()
    for (i in intensity_columns){
      new <- basename(i)
      new_names <- c(new_names, new)
    }

    colnames(df)[length(colnames(df)) - (length(intensity_columns)-1):length(colnames(intensity_columns))] <- new_names
    intensity_names <- grep("\\.(d|raw)$", colnames(df), value = TRUE)
    df[intensity_names] <- sapply(df[intensity_names], as.numeric)
    LOG2.names <- sub("\\.(d|raw)$", ".LOG2", intensity_names)
    df[LOG2.names] <- log2(df[intensity_names])


    ##############################################################################################################
    DIANN_report <- fread(file.path(directory_path, "report.tsv"), sep = "\t", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
    samples <- unique(DIANN_report$Run)[order(unique(DIANN_report$Run))]
    DIANN_report <- DIANN_report %>% #Filtramos por FDR del 1%
      filter(Q.Value < 0.01)

    #Obtenemos los peptidos y peptidos unicos de cada muestra almacenado en listas
    process_subset <- function(df) {
      colnames(df)[5] <- "Unique.Peptide.Count"
      data <- unique(df[, c('Protein.Group', 'Stripped.Sequence')])
      data.all <- df[, c('Protein.Group', 'Stripped.Sequence')]
      unique_peptides_info <- as.data.frame(table(data$Protein.Group))
      peptides_info <- as.data.frame(table(data.all$Protein.Group))

      return(list(data = data, data_all = data.all, unique_peptides_info = unique_peptides_info, peptides_info = peptides_info))
    }

    subsets <- list()
    processed_results <- list()

    #Loop through each unique sample name and create a subset
    for (sample in samples) {
      # Subset the dataframe
      subset_df <- DIANN_report[DIANN_report$Run == sample, ]
      # Store the subset in the list with the sample name as the key
      subsets[[sample]] <- subset_df

      processed_results[[sample]] <- process_subset(subset_df)
    }


    #Extract the dataframes of peptide counts and rename columns accordingly
    peptide_counts_list <- list()
    unique_peptide_counts_list <- list()

    for (sample_name in names(processed_results)) {
      # Extract the unique_peptides_info dataframe
      df_unique_peptides <- processed_results[[sample_name]]$unique_peptides_info #Nos quedamos con el dataframe de cada muestra

      # Rename the columns to indicate the sample name
      colnames(df_unique_peptides) <- c("Var", sample_name)

      # Store the dataframe in the list
      unique_peptide_counts_list[[sample_name]] <- df_unique_peptides
    }

    for (sample_name in names(processed_results)) {
      # Extract the unique_peptides_info dataframe
      df_peptides <- processed_results[[sample_name]]$peptides_info #Nos quedamos con el dataframe de cada muestra

      # Rename the columns to indicate the sample name
      colnames(df_peptides) <- c("Var", sample_name) #Renombramos

      # Store the dataframe in the list
      peptide_counts_list[[sample_name]] <- df_peptides #Guardamos en un dataframe mayor el dataframe de cada muestra
    }
    #Merge all dataframes by the "Var" column
    merged_unique_peptide_counts <- Reduce(function(x, y) merge(x, y, by = "Var", all = TRUE), unique_peptide_counts_list) #Los juntamos todos por la columna Var
    merged_peptide_counts <- Reduce(function(x, y) merge(x, y, by = "Var", all = TRUE), peptide_counts_list) #Los juntamos todos por la columna Var

    #Total_peptides <- unique(DIANN_report[, c('Genes', 'Stripped.Sequence')])
    #peptide_counts_df <- Total_peptides %>%
    #  group_by(Genes) %>%
    #  summarise(Total_Peptides = n(), .groups = "drop")

    #colnames(peptide_counts_df)[1] <- "Var"

    #merged_peptide_counts <- merge(merged_peptide_counts, peptide_counts_df, by = "Var", all = FALSE)

    standardize_columns <- function(df, all_runs) {

      test <- sub("\\.LOG2$", "", all_runs)
      # 1. Replace NAs with 0 (created by the Reduce/Merge step)
      df[is.na(df)] <- 0

      # 2. Identify missing columns
      # We look at column names starting from index 2 (skipping "Var")
      current_runs <- colnames(df)[-1]
      missing_runs <- setdiff(test, current_runs)

      # 3. Add missing columns filled with 0
      if(length(missing_runs) > 0) {
        for(run in missing_runs) {
          df[[run]] <- 0
        }
      }

      return(df)
    }


    # Fix Unique Counts
    merged_unique_peptide_counts <- standardize_columns(merged_unique_peptide_counts, LOG2.names)

    # Fix Peptide Counts
    merged_peptide_counts <- standardize_columns(merged_peptide_counts, LOG2.names)


    colnames(merged_unique_peptide_counts)[2:(length(LOG2.names) +1)] <- paste("Unique peptides", colnames(merged_unique_peptide_counts)[2:(length(LOG2.names) +1)], sep=" ")
    colnames(merged_unique_peptide_counts)[1] <- "Protein"
    colnames(merged_peptide_counts)[2:(length(LOG2.names) +1)] <- paste("Peptides", colnames(merged_peptide_counts)[2:(length(LOG2.names) +1)], sep=" ")
    colnames(merged_peptide_counts)[1] <- "Protein"

    merged_df <- merge(df, merged_unique_peptide_counts, by = "Protein", all = FALSE)

    df <- merge(merged_df, merged_peptide_counts, by = "Protein", all = FALSE)


    # Sample names
    first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]

    condition1_names <- metadata %>%
      filter(group == first_group) %>%
      pull(log2_col)

    second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]


    condition2_names <- metadata %>%
      filter(group == second_group) %>%
      pull(log2_col)


    # Calculate non-NA counts for each condition separately
    df <- df %>%  #Se coge la variable del codigo
      rowwise() %>%
      mutate(
        counts_condition1 = sum(!is.na(c_across(all_of(condition1_names)))),
        counts_condition2 = sum(!is.na(c_across(all_of(condition2_names))))
      ) %>%
      ungroup()

    # Combine the counts into a single non_na_counts column if needed
    df <- df %>%
      mutate(counts = counts_condition1 + counts_condition2)

    # If you need to replace only when non_na_counts is 0 in the original column
    df <- df %>%
      mutate(counts = if_else(counts == 0, counts_condition1 + counts_condition2, counts))
    ##############################################################################################################

    return(df)

  } else if (platform == 4){
    df <- raw

    abundance_names <- grep("Abundance:", intensity_names, value = TRUE)
    df[abundance_names] <- sapply(df[abundance_names], as.numeric)
    LOG2.names <- sub("^Abundance:", "LOG2", abundance_names)
    df[LOG2.names] <- log2(df[abundance_names])

    colnames(df)[colnames(df) == "Accession"] <- "Protein"

    colnames(df)[colnames(df) == "Description"] <- "Protein_description"

    if (organism == 1){
      regex <- regexpr(".*(?=.CGDID)", df$Protein_description, perl = TRUE)
      df$Protein <- regmatches(df$Protein_description, regex)

      #Obtenemos la descripci?n proteica
      regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Protein_description, perl = TRUE)
      df$Protein_description <- regmatches(df$Protein_description, regex3)
    }


  } else if (platform == 5){
    df <- raw

    colnames(df)[1] <- "Protein"
    colnames(df)[4] <- "Protein_description"


    df[intensity_names] <- sapply(df[intensity_names], as.numeric)
    LOG2.names <- sub("Intensity", ".LOG2", intensity_names)
    df[LOG2.names] <- log2(df[intensity_names])

    return(df)

  }

  return(df)


}

##################################################################################
#Obtaining LOG2 names
#' @title Obtain LOG2 Column Names
#' @description A helper function to find all column names containing "LOG2".
#'
#' @param df A data.frame to search.
#'
#' @return A character vector of column names.
#' @export
#'
#' @examples
#' # Create a small example data frame
#' my_df <- data.frame(
#'   Protein = c("P1", "P2"),
#'   Intensity.A = c(100, 200),
#'   LOG2.A = c(6.6, 7.6),
#'   Intensity.B = c(110, 210),
#'   LOG2.B = c(6.7, 7.7)
#' )
#'
#' # This example is runnable!
#' log_names <- obtain_LOG.names(my_df)
#' print(log_names)
#'

obtain_LOG.names <- function(df){

  LOG2.names <- grep("LOG2", colnames(df), value = TRUE)
  return(LOG2.names)
}

##################################################################################
#Subseting unique proteins for each condition
#' @title Obtain Unique Proteins Between Conditions
#' @description Identifies proteins that are uniquely present in one of two
#' conditions, based on the presence of finite (e.g., LOG2 intensity) vs.
#' non-finite (e.g., Inf, NA) values.
#'
#' @details This function is designed for label-free data where a protein's
#' absence is marked by a non-finite value (like Inf from a log-transformation
#' of 0). It identifies proteins present in at least half of the replicates of
#' one condition while being completely absent (non-finite) in the other.
#'
#' @param df A data frame containing expression data, with LOG2 intensity
#'           columns and a 'Protein' column.
#' @param metadata A data frame with sample metadata. Must have 'group' and
#'                 'log2_col' columns.
#' @param selected_conditions A character vector of length 2 specifying the two
#'                            groups from the 'group' column to compare.
#'
#' @return A list containing two data frames:
#'         1. `cond1_unicas`: A data frame of proteins unique to the first group.
#'         2. `cond2_unicas`: A data frame of proteins unique to the second group.
#' @export
#' @importFrom dplyr %>% filter pull
#'
#' @examples
#' # Create synthetic data with infinite values (representing missing)
#' df_prot <- data.frame(
#'   Protein = paste0("Prot", 1:4),
#'   LOG2.C1 = c(10, 12, Inf, 9),
#'   LOG2.C2 = c(11, 13, Inf, 10),
#'   LOG2.T1 = c(Inf, Inf, 8, 11),
#'   LOG2.T2 = c(Inf, Inf, 7, 12)
#' )
#' # Prot1/2 are present in Control, absent in Treatment
#' # Prot3 is present in Treatment, absent in Control
#' # Prot4 is present in both
#'
#' metadata <- data.frame(
#'   log2_col = c("LOG2.C1", "LOG2.C2", "LOG2.T1", "LOG2.T2"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' # Note: T is [1] (second_group), C is [2] (first_group)
#' selected_cond <- c("Treatment", "Control")
#'
#' unique_list <- obtain_unique_proteins(df_prot, metadata, selected_cond)
#'
#' # cond1_unicas (Control) should have Prot1, Prot2
#' print("Unique to Control:")
#' print(unique_list[[1]]$Protein)
#'
#' # cond2_unicas (Treatment) should have Prot3
#' print("Unique to Treatment:")
#' print(unique_list[[2]]$Protein)
#'

obtain_unique_proteins <- function(df, metadata, selected_conditions){

  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]

  condition1_names <- metadata %>%
    dplyr::filter(group == first_group) %>%
    dplyr::pull(log2_col)

  replicas_condicion1 <- length(condition1_names)
  print(replicas_condicion1)

  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]

  condition2_names <- metadata %>%
    dplyr::filter(group == second_group) %>%
    dplyr::pull(log2_col)

  replicas_condicion2 <- length(condition2_names)
  print(replicas_condicion2)


  condi_names <- c(condition1_names,condition2_names)

  df2 <- df[condi_names]
  df2 <- as.matrix(df2)

  total_replicas <- replicas_condicion1 + replicas_condicion2

  finite_sums <- rowSums(is.finite(df2[, 1:replicas_condicion1]))
  infinite_sums <- rowSums(!is.finite(df2[, (replicas_condicion1 + 1):total_replicas]))

  cond1_exclusive <- Reduce(`|`, lapply(0:(replicas_condicion2 - 1), function(i) {
    condition <- finite_sums == (replicas_condicion1 - i) & infinite_sums == replicas_condicion2
    sufficient_presence <- finite_sums >= ceiling(replicas_condicion1 / 2)
    condition & sufficient_presence
  }))

  cond2_exclusive <- Reduce(`|`, lapply(1:replicas_condicion2, function(i) {
    condition <- finite_sums == 0 & infinite_sums == (replicas_condicion2 - i)
    sufficient_presence <- infinite_sums <= ceiling(replicas_condicion2 / 2)
    condition & sufficient_presence
  }))

  df$cond1_exclusive <- cond1_exclusive
  df$cond2_exclusive <- cond2_exclusive

  cond1_unicas <- dplyr::filter(df, df$cond1_exclusive)
  cond2_unicas <- dplyr::filter(df, df$cond2_exclusive)

  return(list(cond1_unicas, cond2_unicas))
}

##################################################################################
#Proteins identified


#' @title Identify Proteins
#' @description Identifies unique and common proteins between two data frames.
#'
#' @param df1 A data frame for condition 1, must have a "Protein" column.
#' @param df2 A data frame for condition 2, must have a "Protein" column.
#'
#' @return A list containing three data frames: `unique_df1` (proteins unique
#' to df1), `unique_df2` (proteins unique to df2), and `common_df` (proteins
#' common to both).
#' @export
#' @importFrom dplyr bind_rows distinct select
#'
#' @examples
#' df1 <- data.frame(Protein = c("P1", "P2", "P3"), Value = 1:3)
#' df2 <- data.frame(Protein = c("P2", "P3", "P4"), Value = 2:4)
#' id_list <- identify_proteins(df1, df2)
#' print(id_list$common_df)
#' print(id_list$unique_df1)
#'
identify_proteins <- function(raw, metadata, platform, selected_conditions) {

  # --- 1. Get sample names for each condition ---
  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]
  condition1_names <- metadata %>%
    dplyr::filter(group == first_group) %>%
    dplyr::pull(intensity_sample_name)

  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]
  condition2_names <- metadata %>%
    dplyr::filter(group == second_group) %>%
    dplyr::pull(intensity_sample_name)

  condi.names <- c(condition1_names, condition2_names)

  # --- 2. Filter data based on platform ---
  filt <- list()
  if (platform %in% c(1, 2, 4, 5)) {
    for (i in 1:length(condi.names)) {
      filt[[i]] <- raw[raw[, condi.names[i]] > 0, ]
    }
  } else if (platform == 3) {
    for (i in 1:length(condi.names)) {
      filt[[i]] <- raw[complete.cases(raw[, condi.names[i]]), ]
    }
  }

  protein_counts <- sapply(filt, nrow)
  plot_data <- data.frame(
    Sample = condi.names,
    Count = protein_counts,
    Condition = rep(c(first_group, second_group), times = c(length(condition1_names), length(condition2_names)))
  )

  plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Count, fill = Condition)) +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::scale_fill_manual(values = c("light green", "light blue")) + # You can customize colors here
    ggplot2::labs(
      title = "Proteínas Cuantificadas",
      x = "Muestras",
      y = "Total de Proteínas"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1), # Better angle for long names
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(plot_object)
}
#Venn Diagram
#' @title Create a Venn Diagram
#' @description Generates a Venn diagram plot object from a list of unique
#'              proteins and a data frame of common proteins.
#'
#' @param df A data frame of common proteins (must have a "Protein" column).
#' @param unique_proteins A list containing two data frames:
#'                        1. Proteins unique to the first set.
#'                        2. Proteins unique to the second set.
#' @param label1 Label for the first set (e.g., "Control").
#' @param label2 Label for the second set (e.g., "Treatment").
#' @param color1 Color for the first set.
#' @param color2 Color for the second set.
#'
#' @return A gList object (a plot) representing the Venn diagram.
#'
#' @importFrom VennDiagram venn.diagram
#' @importFrom dplyr bind_rows distinct
#' @import grid
#' @export
#'
#' @examples
#' # --- Create 100% Synthetic Data ---
#' # 1. Create a data frame of proteins common to both groups
#' common_proteins <- data.frame(
#'   Protein = c("ProteinA", "ProteinB", "ProteinC")
#' )
#'
#' # 2. Create a list of data frames for unique proteins
#' unique_protein_list <- list(
#'   # Proteins unique to Group 1 ("Control")
#'   data.frame(Protein = c("ProteinD", "ProteinE")),
#'   # Proteins unique to Group 2 ("Treatment")
#'   data.frame(Protein = c("ProteinF", "ProteinG", "ProteinH", "ProteinI"))
#' )
#'
#' # --- Run the Example ---
#' # This code is now self-contained and runnable by BiocCheck
#' my_venn_plot <- venn_diagram(
#'   df = common_proteins,
#'   unique_proteins = unique_protein_list,
#'   label1 = "Control",
#'   label2 = "Treatment",
#'   color1 = "blue",
#'   color2 = "red"
#' )
#'
#' # Venn diagrams often need to be explicitly drawn
#' grid.newpage()
#' grid.draw(my_venn_plot)
#'
venn_diagram <- function(df, unique_proteins, label1, label2, color1, color2){

  # Create cond1_set by adding unique_proteins[[2]] to df
  cond1_set <- dplyr::bind_rows(df, unique_proteins[[2]]) %>% dplyr::distinct()

  # Create cond2_set by adding unique_proteins[[1]] to df
  cond2_set <- dplyr::bind_rows(df, unique_proteins[[1]]) %>% dplyr::distinct()

  # Create sets for Venn diagram using Proteins
  x <- list(cond1_set$Protein, cond2_set$Protein)
  names(x) <- c(label1, label2)

  # Generate Venn diagram
  venn.plot <- VennDiagram::venn.diagram(
    x = x,
    category.names = c(label1, label2),
    filename = NULL,
    output = FALSE,
    imagetype = "png",
    main = "Venn Diagram",
    fill = c(color1, color2),
    disable.logging = TRUE
  )

  return(venn.plot)
}


##################################################################################
#Filtering
#' @title Filter proteins by valid values
#' @description A helper function to filter a protein expression data frame
#' (e.g., from MaxQuant) based on a minimum number of valid
#' (non-NA, non-zero) values across samples.
#'
#' @param raw_data A data frame containing the raw protein expression data.
#' @param metadata A data frame mapping intensity column names to groups.
#' @param valid_thr A numeric threshold for the minimum number of valid
#'                  values. For example, 2 means a protein must be valid
#'                  in at least 2 samples to be kept.
#'
#' @return A list containing the filtered data frame ("filtered") and
#'         the LOG2-transformed data frame ("LOG2").
#' @export
#'
#' @examples
#' # Create TINY, SYNTHETIC data for the example:
#' my_raw_data <- data.frame(
#'    Protein.ID = c("P1", "P2", "P3"),
#'    Intensity.S1 = c(100, 200, 0),
#'    Intensity.S2 = c(110, 210, 300),
#'    Potential.contaminant = c("", "", ""),
#'    Reverse = c("", "", ""),
#'    Only.identified.by.site = c("", "", "")
#' )
#'
#' my_metadata <- data.frame(
#'    intensity_sample_name = c("Intensity.S1", "Intensity.S2"),
#'    group = c("A", "B"),
#'    log2_col = c("LOG2.S1", "LOG2.S2")
#' )
#'
#' # This example is self-contained and runnable
#' filtered_list <- filter_valids(my_raw_data, my_metadata, 2)
#' print(filtered_list$filtered)
#'
filter_valids <- function(df, metadata, unique_proteins, min_prop = NULL,
                          at_least_one = FALSE, labeltype = 1) {

  groups <- unique(metadata$group)
  if (length(groups) < 2) stop("You need at least two conditions.")

  cond.names <- lapply(groups, function(g) {
    metadata %>%
      dplyr::filter(group == g) %>%
      dplyr::pull(log2_col)
  })
  names(cond.names) <- groups

  min_count <- sapply(cond.names, function(cols) ceiling(length(cols) * min_prop))

  # Remove unique proteins
  total_unique_proteins <- rbind(unique_proteins[[1]], unique_proteins[[2]])
  if (nrow(total_unique_proteins) > 0) {
    common_df <- df %>%
      dplyr::filter(!Protein %in% total_unique_proteins$Protein)
  } else {
    common_df <- df
  }

  if (labeltype == 1) {
    # Filter by non-NA counts per condition
    cond.filter <- sapply(seq_along(cond.names), function(i) {
      mat <- as.matrix(common_df[cond.names[[i]]])
      rowSums(is.finite(mat)) >= min_count[i]
    })

    common_df$KEEP <- if (at_least_one) {
      apply(cond.filter, 1, any)
    } else {
      apply(cond.filter, 1, all)
    }

    common_df <- common_df %>% dplyr::filter(KEEP)

    # Replace non-finite values with NA
    cols_to_modify <- unlist(cond.names)
    common_df[cols_to_modify] <- lapply(common_df[cols_to_modify], function(col) {
      col[!is.finite(col)] <- NA
      return(col)
    })
  } else if (labeltype == 2) {
    # Label-based: remove missing (assumes external function)
    common_df <- df
    common_df <- remove_missing(common_df) # Assuming remove_missing is your own function
  }

  # ----------- Vectorized CV Calculation (for both label types) -------------
  condition1_names <- cond.names[[1]]
  condition2_names <- cond.names[[2]]

  valuesA <- as.matrix(common_df[, condition1_names])
  valuesB <- as.matrix(common_df[, condition2_names])
  valuesC <- cbind(valuesA, valuesB)

  row_cv <- function(mat) {
    row_means <- rowMeans(mat, na.rm = TRUE)
    # 'sd' is from base R (stats package)
    row_sds <- apply(mat, 1, sd, na.rm = TRUE)
    cv <- (row_sds / row_means) * 100
    cv[is.nan(cv)] <- NA
    return(cv)
  }

  common_df$CV_Control <- row_cv(valuesA)
  common_df$CV_Tratamiento <- row_cv(valuesB)
  common_df$CV_TratamientoyControl <- row_cv(valuesC)
  # ---------------------------------------------------------------------------

  return(common_df)
}


#' @title Filter Proteins by Unique Peptide Counts
#' @description Filters a protein data frame based on a minimum number of
#' unique peptides found in a minimum fraction of replicates per condition.
#'
#' @details This function is useful for filtering data (e.g., from MaxQuant) to
#' remove proteins with low confidence. It keeps a protein if, for
#' **at least one** experimental group (condition), the unique peptide count is
#' `>= number` in at least `min_fraction` of the samples in that group.
#'
#' @param df A data frame with protein data. Must have rownames
#'           (e.g., Protein IDs) and columns corresponding to unique
#'           peptide counts.
#' @param metadata A metadata data frame. Must contain 'group' and
#'                 'unique_peptides_col' columns to map samples to groups.
#' @param number Numeric. The minimum number of unique peptides (e.g., 1 or 2)
#'               required to count as 'present' in a single sample.
#' @param min_fraction Numeric. The minimum fraction (0 to 1) of replicates
#'                     in *at least one* group that must meet the 'number'
#'                     threshold. (e.g., 0.5 means at least 50% of replicates).
#'
#' @return A data frame (`df_filtered`) containing only the rows (proteins)
#'         that pass the filter.
#' @export
#' @importFrom dplyr %>% filter pull
#'
#' @examples
#' # Create synthetic data
#' peptide_df <- data.frame(
#'   Pep.C1 = c(2, 0, 1, 5),
#'   Pep.C2 = c(3, 1, 0, 6),
#'   Pep.T1 = c(0, 2, 1, 0),
#'   Pep.T2 = c(0, 3, 0, 0)
#' )
#' rownames(peptide_df) <- c("ProtA", "ProtB", "ProtC", "ProtD")
#'
#' meta <- data.frame(
#'   unique_peptides_col = c("Pep.C1", "Pep.C2", "Pep.T1", "Pep.T2"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' # Filter: require min 1 peptide in at least 75% of reps (i.e., 2/2)
#' # ProtA: Passes in Control (2/2 have >= 1). Kept.
#' # ProtB: Passes in Treatment (2/2 have >= 1). Kept.
#' # ProtC: Fails in Control (1/2), Fails in Treatment (1/2). Dropped.
#' # ProtD: Passes in Control (2/2 have >= 1). Kept.
#'
#' filtered_df <- unique_peptides_filter(
#'   df = peptide_df,
#'   metadata = meta,
#'   number = 1,
#'   min_fraction = 0.75
#' )
#'
#' print("Filtered proteins (should be A, B, D):")
#' print(rownames(filtered_df))
#'
unique_peptides_filter <- function(df, metadata, number, min_fraction) {

  # Get sample column names per condition
  groups <- unique(metadata$group)
  if (length(groups) < 2) stop("You need at least two conditions.")

  condition_columns <- lapply(groups, function(g) {
    metadata %>%
      dplyr::filter(group == g) %>%
      dplyr::pull(unique_peptides_col)
  })
  names(condition_columns) <- groups

  # Function to filter rows (proteins) based on min_fraction
  filter_proteins_by_condition <- function(data_subset, condition_name) {
    num_samples <- ncol(data_subset)
    min_required <- ceiling(min_fraction * num_samples)
    rowSums(data_subset >= number, na.rm = TRUE) >= min_required
  }

  # Apply filtering per condition
  filtered_proteins_list <- mapply(function(cols, cond_name) {
    subset <- df[, cols, drop = FALSE]
    keep_rows <- filter_proteins_by_condition(subset, cond_name)
    df[keep_rows, , drop = FALSE]
  }, condition_columns, names(condition_columns), SIMPLIFY = FALSE)

  # Get union of proteins passing any condition
  all_filtered_proteins <- Reduce(union, lapply(filtered_proteins_list, rownames))

  if (length(all_filtered_proteins) == 0) {
    warning("No proteins meet the filtering criteria.")
    return(data.frame())
  }

  df_filtered <- df[rownames(df) %in% all_filtered_proteins, ]
  return(df_filtered)
}



##################################################################################
#Normalization
#' @title Perform median centering normalization
#' @description Normalizes a numeric matrix of log-expression values by
#' subtracting the median of each column.
#'
#' @param df A numeric matrix where rows are proteins and columns are samples.
#'
#' @return A numeric matrix of the same dimensions, now median-centered.
#' @export
#'
#' @examples
#' # Create a small, synthetic log-expression matrix
#' log_matrix <- matrix(c(10, 10.2, 9.8,  # Sample 1
#'                        12, 12.1, 11.9), # Sample 2
#'                      nrow = 3, ncol = 2,
#'                      dimnames = list(c("ProtA", "ProtB", "ProtC"),
#'                                      c("Sample1", "Sample2")))
#'
#' print("Original Matrix:")
#' print(log_matrix)
#'
#' # Run the normalization
#' normalized_matrix <- median_centering(log_matrix)
#'
#' print("Normalized Matrix:")
#' print(normalized_matrix)
#'
#' # The new medians should be 0
#' print("New Column Medians (should be 0):")
#' print(apply(normalized_matrix, 2, median))
#'
median_centering <- function(df, LOG2.names) {

  df[, LOG2.names] <- lapply(LOG2.names,
                             function(x) {
                               LOG2 <- df[[x]]
                               LOG2[!is.finite(LOG2)] <- NA
                               gMedian <- median(LOG2, na.rm = TRUE)
                               LOG2 - gMedian #Se calcula el valor.
                             }
  )

  return(df)
}
#' @title Apply Normalization to LOG2 Columns
#' @description A wrapper function to apply various normalization methods from
#' the `wrMisc` package to the specified LOG2 intensity columns.
#'
#' @param df A data frame containing the LOG2 expression data.
#' @param LOG2.names A character vector of column names (from `df`) to which
#'                   the normalization should be applied.
#' @param method A character string specifying the normalization method to be
#'               passed to `wrMisc::normalizeThis`. Examples include
#'               "median", "mean", "trimMean", "quant", etc.
#'
#' @return The original data frame (`df`) with the specified `LOG2.names`
#'         columns replaced by their normalized versions.
#' @export
#' @importFrom wrMisc normalizeThis
#'
#' @examples
#' # Create synthetic LOG2 data
#' log_df <- data.frame(
#'   Protein = c("ProtA", "ProtB"),
#'   LOG2.Sample1 = c(10, 12),
#'   LOG2.Sample2 = c(11, 13)
#' )
#'
#' # Specify which columns to normalize
#' log_cols <- c("LOG2.Sample1", "LOG2.Sample2")
#'
#' # Run median normalization
#' normalized_df <- normalization_func(log_df, log_cols, "median")
#'
#' print("Original Data:")
#' print(log_df)
#' print("Normalized Data:")
#' print(normalized_df)
#'
normalization_func <- function(df, LOG2.names, method){
  if (method == "trimMean"){
    df[, LOG2.names] <- wrMisc::normalizeThis(dat = df[, LOG2.names], method = method, trimFa = 0.4)
  } else if (method != "trimMean"){
    df[, LOG2.names] <- wrMisc::normalizeThis(dat = df[, LOG2.names], method = method)
  }
  return(df)
}



##################################################################################
#Imputation
#' @title Impute Missing Values using Down-Shifted Normal Distribution
#' @description Replaces non-finite (NA, Inf) values in specified LOG2 columns
#' with random values drawn from a normal distribution. This distribution is
#' shifted to a lower mean and narrower standard deviation relative to the
#' observed valid data.
#'
#' @details This method is a common "down-shift" or "Perseus-style" imputation
#' for label-free proteomics, simulating the low-abundance signals of
#' proteins near the detection limit. It assumes missingness is not at random
#' (MNAR).
#'
#' @param df A data frame. This frame **must** contain a logical column named
#'           `KEEP` which specifies the rows to be used for calculating the
#'           mean and standard deviation of the valid data.
#' @param LOG2.names A character vector of column names (from `df`) to impute.
#' @param width Numeric. The factor by which to shrink the standard deviation
#'              of the valid data (e.g., 0.3).
#' @param downshift Numeric. The factor (in standard deviations) by which to
#'                  shift the mean of the valid data (e.g., 1.8).
#'
#' @return The data frame `df` with imputed LOG2 columns and new 'impute.*'
#'         logical columns indicating which values were imputed.
#' @export
#' @importFrom stats sd rnorm
#'
#' @examples
#' # Create synthetic LOG2 data with non-finite values (Inf)
#' log_df <- data.frame(
#'   Protein = c("ProtA", "ProtB", "ProtC", "ProtD"),
#'   LOG2.S1 = c(10, 12, Inf, 8),
#'   LOG2.S2 = c(11, 13, Inf, Inf),
#'   # Add a KEEP column (e.g., from a filter)
#'   # Here, we keep all proteins for calculation
#'   KEEP = c(TRUE, TRUE, TRUE, TRUE)
#' )
#'
#' # Specify columns to impute
#' log_cols <- c("LOG2.S1", "LOG2.S2")
#'
#' # Run imputation
#' set.seed(123) # for reproducible example
#' imputed_df <- impute_data(log_df, log_cols)
#'
#' print("Original Data:")
#' print(log_df)
#' print("Imputed Data:")
#' print(imputed_df)
#'
impute_data <- function(df, LOG2.names, width = 0.3, downshift = 1.8) {

  impute.names <- sub("LOG2", "impute", LOG2.names)

  # Create new columns indicating whether the values are imputed
  df[impute.names] <- lapply(LOG2.names, function(x) !is.finite(df[, x]))

  # Imputation

  df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- df[[x]]
                             temp[!is.finite(temp)] = NA

                             temp.sd <- width * sd(temp[df$KEEP], na.rm = TRUE) # shrink sd width
                             temp.mean <- mean(temp[df$KEEP], na.rm = TRUE) -
                               downshift * sd(temp[df$KEEP], na.rm = TRUE) # shift mean of imputed values

                             n.missing <- sum(is.na(temp))
                             temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
                             return(temp)
                           })


  return(df)
}


#' @title Impute missing data using k-NN
#' @description Uses K-Nearest Neighbors (k-NN) to impute missing (NA)
#' values in a log-expression matrix.
#'
#' @param df A numeric matrix with NA values.
#'
#' @return A numeric matrix with NA values imputed.
#' @importFrom VIM kNN
#' @export
#'
#' @examples
#' # Create a matrix with missing (NA) values
#' log_matrix_na <- matrix(c(10, 10.2, NA,   # Sample 1 (has NA)
#'                           12, NA, 11.9),  # Sample 2 (has NA)
#'                         nrow = 3, ncol = 2,
#'                         dimnames = list(c("ProtA", "ProtB", "ProtC"),
#'                                         c("Sample1", "Sample2")))
#'
#' print("Matrix with NAs:")
#' print(log_matrix_na)
#'
#' # Run the imputation
#' # This requires the VIM package, which we import
#' imputed_matrix <- impute_knn_data(log_matrix_na)
#'
#' print("Imputed Matrix (NAs removed):")
#' print(imputed_matrix)


impute_KNN_data <- function(df, LOG2.names, ...){

  impute.names <- sub("LOG2", "impute", LOG2.names)
  df[impute.names] <- lapply(LOG2.names, function(x) !is.finite(df[, x]))

  #Imputaci?n
  df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- df[[x]]
                             temp[!is.finite(temp)] = NA
                             return(temp)

                           })

  imp_knn <- VIM::kNN(df, variable = LOG2.names, ... )

  return(imp_knn)
}




##################################################################################
#Processing check

#' @title Plot Coefficient of Variation vs. Mean
#' @description Generates a plot of the squared Coefficient of Variation (CV)
#' versus the mean abundance (A). This is often used to check imputation
#' quality.
#'
#' @param y A numeric data frame or matrix where rows are features (proteins)
#'          and columns are samples (LOG2 intensities).
#' @param trend Logical. If `TRUE`, a `limma::loessFit` trendline is added
#'              to the plot.
#' @param main Character. The main title for the plot.
#' @param ... Additional arguments passed to the `plot` function.
#'
#' @return A data frame containing the `mean` (A) and `CV` (CV^2) values
#'         that were plotted.
#' @export
#' @importFrom stats na.omit
#' @importFrom graphics plot lines
#' @importFrom matrixStats rowSds
#' @importFrom limma loessFit
#'
#' @examples
#' # Create synthetic LOG2 intensity data
#' log_data <- matrix(rnorm(100, mean = 10, sd = 1), nrow = 20, ncol = 5)
#' log_data[sample(1:100, 10)] <- NA # Add some NAs
#'
#' # Create the CV plot
#' cv_data <- plotCV2(log_data, main = "CV vs. Mean Plot")
#' print(head(cv_data))
#'
#'
plotCV2 <- function(y, trend = TRUE, main= "Imputation check", ...){
  y <- na.omit(y) # na.omit is base R (stats)
  A <- rowMeans(y, na.rm = TRUE) # rowMeans is base R
  CV <- (matrixStats::rowSds(data.matrix(y), na.rm = TRUE)/A)^2 # data.matrix is base R
  res <- data.frame(mean = A, CV = CV)
  plot(A, CV ,  ylim = c(min(CV)-0.001, max(CV) +0.001),  ...) # plot, min, max are base R
  if(trend){
    fit <- limma::loessFit(CV, A)
    o <- order(A) # order is base R
    lines(A[o], fit$fitted[o], lwd =2, col = "red") # lines is base R
  }

  return(res)
}

#' @title Create Boxplot of Sample Intensities
#' @description Generates a ggplot2 boxplot of LOG2 intensities, grouped
#'              by sample and colored by condition.
#'
#' @param df A data frame containing the LOG2 intensity data.
#' @param metadata A metadata data frame. Must contain 'group' and
#'                 'log2_col' columns to map samples to groups.
#' @param selected_conditions A character vector of length 2 specifying the
#'                            two conditions to compare.
#'
#' @return A ggplot object.
#' @export
#' @importFrom dplyr %>% filter pull select all_of everything mutate case_when
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_classic theme
#' @importFrom ggplot2 element_text
#' @importFrom stats na.omit
#'
#' @examples
#' # Create synthetic data
#' log_df <- data.frame(
#'   LOG2.C1 = rnorm(20, 10),
#'   LOG2.C2 = rnorm(20, 10.1),
#'   LOG2.T1 = rnorm(20, 12),
#'   LOG2.T2 = rnorm(20, 12.1)
#' )
#'
#' meta <- data.frame(
#'   group = c("Control", "Control", "Treatment", "Treatment"),
#'   log2_col = c("LOG2.C1", "LOG2.C2", "LOG2.T1", "LOG2.T2")
#' )
#'
#' # Generate the boxplot
#' boxplot_function(log_df, meta, selected_conditions = c("Control", "Treatment"))
#'
boxplot_function <- function(df, metadata, selected_conditions) {

  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]
  condition1_names <- metadata %>%
    dplyr::filter(group == first_group) %>%
    dplyr::pull(log2_col)

  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]
  condition2_names <- metadata %>%
    dplyr::filter(group == second_group) %>%
    dplyr::pull(log2_col)

  cond.names <- c(condition1_names, condition2_names)

  plot_data_long <- df %>%
    dplyr::select(dplyr::all_of(cond.names)) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    dplyr::mutate(Condition = dplyr::case_when(
      Sample %in% condition1_names ~ first_group,
      Sample %in% condition2_names ~ second_group
    )) %>%
    na.omit() # Remove rows with missing intensities (na.omit is base R)

  plot_object <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = Sample, y = Intensity, fill = Condition)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      title = "Boxplot of Intensities",
      y = "Log2(Intensity)",
      x = "Sample"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(plot_object)
}


#' @title Plot Missingness Pattern
#' @description Generates a missing value aggregation plot from the `VIM`
#' package to visualize patterns of missing data.
#'
#' @param df A data frame containing the data to be plotted.
#' @param cond.names A character vector of column names from `df` to be
#'                   included in the missingness plot.
#'
#' @return A plot showing the aggregation of missing values.
#' @export
#' @importFrom graphics par
#' @importFrom VIM aggr
#'
#' @examples
#' # Create synthetic data with NAs
#' log_df <- data.frame(
#'   LOG2.C1 = c(10, 10.1, NA, 9.9),
#'   LOG2.C2 = c(11, NA, NA, 11.2),
#'   LOG2.T1 = c(12, 12.1, 12.2, 12.3)
#' )
#'
#' # Plot the missingness pattern
#' preimputation_state(log_df, c("LOG2.C1", "LOG2.C2", "LOG2.T1"))
#'
#'
preimputation_state <- function(df, cond.names){

  par(mfrow=c(2,2)) # par is base R
  VIM::aggr(df[, cond.names], delimiter = "NA", labels = names(df[cond.names]), cex.axis = 0.4)
}


#' @title Plot Post-Imputation Density
#' @description Plots the density of intensity values before and after
#' imputation to visualize the effect of the imputation.
#'
#' @param data_filtered The data frame *before* imputation, containing
#'                      non-finite values (NA, Inf).
#' @param imputation An indicator (1 or 2) specifying the imputation method.
#'                   1 = `impute_KNN_data`, 2 = `impute_data`.
#' @param LOG2.names A character vector of the LOG2 column names to use.
#'
#' @return A ggplot object showing the two density plots.
#' @export
#' @importFrom dplyr %>% mutate bind_rows starts_with
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_density theme_minimal labs
#'
#' @examples
#' # This function depends on `impute_KNN_data` and `impute_data`
#' # We must create those functions and the data they need.
#'
#' # 1. Create synthetic data
#' data_filt <- data.frame(
#'   Protein = c("A", "B", "C", "D"),
#'   LOG2.S1 = c(10, 12, NA, 8),
#'   LOG2.S2 = c(11, 13, 9, NA),
#'   KEEP = c(TRUE, TRUE, TRUE, TRUE) # For impute_data
#' )
#'
#' # 2. Create minimal versions of the imputation functions for the example
#' impute_KNN_data <- function(df, lnames, k) {
#'    df[lnames][is.na(df[lnames])] <- mean(unlist(df[lnames]), na.rm = TRUE)
#'    return(df)
#' }
#' impute_data <- function(df, lnames) {
#'    df[lnames][is.na(df[lnames])] <- rnorm(sum(is.na(df[lnames])), mean = 7, sd = 0.3)
#'    return(df)
#' }
#'
#' # 3. Run the plot function
#' postimputation_state(data_filt, imputation = 2, LOG2.names = c("LOG2.S1", "LOG2.S2"))
#'
postimputation_state <- function(data_filtered, imputation, LOG2.names){

  df_no_imputed <- data_filtered[LOG2.names] #quedarnos con columnas con intensidad LOG2

  if (imputation == 1) {
    data <- impute_KNN_data(as.data.frame(data_filtered), LOG2.names, k = 5) # as.data.frame is base R
  } else if (imputation == 2){
    data <- impute_data(as.data.frame(data_filtered), LOG2.names)

  }

  df_imputed <- data[LOG2.names] #Imputar y lo mismo

  # Convert to long format
  df_imputed_long <- df_imputed %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("LOG2"), names_to = "Sample", values_to = "Intensity") %>%
    dplyr::mutate(Type = "Imputed")

  df_non_imputed_long <- df_no_imputed %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("LOG2"), names_to = "Sample", values_to = "Intensity") %>%
    dplyr::mutate(Type = "Original")


  df_combined <- dplyr::bind_rows(df_non_imputed_long, df_imputed_long)


  ggplot2::ggplot(df_combined, ggplot2::aes(x = Intensity, fill = Type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Intensity Distributions: Imputed vs Original", x = "Intensity", y = "Density")
}


#' @title Plot a Histogram
#' @description A simple wrapper to plot a histogram for a specific column
#'              in a data frame.
#'
#' @param df A data frame.
#' @param colname The name of the (numeric) column to plot.
#' @param color The color for the histogram.
#' @param title The main title for the plot.
#'
#' @return A histogram plot.
#' @export
#' @importFrom graphics hist
#'
#' @examples
#' my_data <- data.frame(p.value = rnorm(100))
#' histogram(my_data, "p.value", "lightblue", "P-Value Distribution")
#'
#'
histogram <- function(df, colname, color, title){
  hist(as.numeric(df[,colname]), main = title, xlab = 'Intensity-value', col = color)


}
#' @title Create a Scatterplot with Marginal Densities
#' @description Uses `ggplot2` and `ggExtra` to create a scatterplot comparing
#'              two columns, with marginal density plots.
#'
#' @param df A data frame.
#' @param colname1 Character. The name of the column for the x-axis.
#' @param colname2 Character. The name of the column for the y-axis.
#'
#' @return A `ggExtra` plot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs
#' @importFrom ggExtra ggMarginal
#'
#' @examples
#' my_data <- data.frame(
#'   Sample1 = rnorm(100),
#'   Sample2 = rnorm(100)
#' )
#' scatterplot_function(my_data, "Sample1", "Sample2")
#'
scatterplot_function <- function(df, colname1, colname2){
  corr_plot <- ggplot2::ggplot(df, ggplot2::aes(df[,colname1], df[,colname2])) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = colname2, x = colname1)

  ggExtra::ggMarginal(corr_plot, type = "densigram")

}


#' @title Create a Q-Q Plot
#' @description Generates a Q-Q plot to compare the distributions of two
#'              columns in a data frame.
#'
#' @param df A data frame.
#' @param colname1 Character. The name of the column for the x-axis.
#' @param colname2 Character. The name of the column for the y-axis.
#' @param color The color for the points.
#'
#' @return A ggplot object.
#' @export
#' @importFrom stats qqnorm
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal geom_qq_line
#' @importFrom ggplot2 ylab xlab ggtitle
#'
#' @examples
#' my_data <- data.frame(
#'   Theoretical = rnorm(100),
#'   Sample = rnorm(100, mean = 0.5)
#' )
#' qqplot_function(my_data, "Theoretical", "Sample", "blue")
#'
qqplot_function <- function(df, colname1, colname2, color){

  p <- data.frame(qqnorm(x = df[,colname1], y = df[,colname2])) # qqnorm is base R (stats)
  ggplot2::ggplot(p, ggplot2::aes(x, y)) +
    ggplot2::geom_point(pch = 21,
                        colour = color) +
    ggplot2::theme_minimal() +
    ggplot2::geom_qq_line(ggplot2::aes(sample = y)) +
    ggplot2::ylab("Sample dist") +
    ggplot2::xlab("Theoretical dist") +
    ggplot2::ggtitle("Q-Q plot")
}


#' Generate a ggplot2 Correlation Matrix Heatmap
#' @title Correlation Matrix
#' @description Calculates a Pearson correlation matrix from a data frame and visualizes it
#' as a tile-based heatmap using `ggplot2`, overlaying the correlation coefficients.
#'
#' @param df A data frame. The columns for which correlations are desired should
#'   be numeric. Non-numeric columns will be ignored by `cor()`.
#' @param display A string specifying the visual representation (e.g., "circle").
#'   **Note:** This parameter is currently not implemented in the function body,
#'   which defaults to 'tile'.
#' @param tl.col A string for text label color (i.e., the variable names).
#'   **Note:** This parameter is currently not implemented in the function body.
#' @param addCoef.col A string specifying the color for the overlaid
#'   correlation coefficients (e.g., `"black"`).
#'
#' @return A `ggplot` object representing the correlation heatmap.
#'
#' @details
#' The function first calculates the correlation matrix using `cor()` with
#' `use = "complete.obs"`, which handles missing values by pairwise deletion.
#' It then melts the matrix into a "long" format suitable for `ggplot2`.
#' The plot uses `scale_fill_gradient2` to create a divergent color scale
#' (Red -> White -> Blue) centered at 0. `coord_fixed()` ensures square tiles.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal
#' @importFrom ggplot2 theme element_text element_blank coord_fixed geom_text
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom stats cor
#'
#' @export
#'
#' @examples
#' # Create a sample data frame
#' df_sample <- data.frame(
#'   A = rnorm(20),
#'   B = rnorm(20),
#'   C = rnorm(20)
#' )
#'
#' # Generate the plot
#' p <- corrplot_function(df_sample, addCoef.col = "blue")
#'
#' # To display the plot, run:
#' # print(p)
#'
corrplot_function <- function(df, display = "circle", tl.col = "black", addCoef.col = "black") {

  df_cor <- cor(df, use = "complete.obs", method = "pearson") # cor is base R

  # Modern way using tidyr and tibble (both part of tidyverse)
  melted_cormat <- as.data.frame(df_cor) %>%
    tibble::rownames_to_column("Var1") %>%
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "value")

  plot_object <- ggplot2::ggplot(data = melted_cormat, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = "red", high = "blue", mid = "white",
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = "Pearson\nCorrelation"
    ) +
    ggplot2::theme_minimal() + # A clean theme
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed() # Ensures the plot is square

  plot_object <- plot_object +
    ggplot2::geom_text(
      ggplot2::aes(label = round(value, 2)), # Round coefficients to 2 decimal places
      color = addCoef.col,
      size = 3
    )
  return(plot_object)
}

##################################################################################
#Differential expression analysis
##################################################################################

#' Perform Differential Expression Statistical Analysis
#' @title Limma function for statistical analysis
#' @description A comprehensive function to perform differential expression (DE) analysis
#' between two conditions. It supports t-tests, Wilcoxon tests, and `limma`,
#' with options for paired data and `DEqMS` spectral count weighting.
#'
#' @param df The primary data frame containing expression data (rows=proteins,
#'   cols=samples). Must also contain 'Protein' and 'Protein_description' columns.
#' @param test Numeric. The statistical test to use: `1` for `t.test`, `2` for
#'   `limma`, `3` for `wilcox.test`.
#' @param paired Logical. `TRUE` if the experimental design is paired, `FALSE`
#'   otherwise.
#' @param metadata A data frame with sample metadata. Must contain `group`
#'   (condition) and `log2_col` (sample names matching `df` colnames).
#' @param logfcup Numeric. The log-fold change (LFC) threshold for
#'   'Up-regulated' status.
#' @param logfcdown Numeric. The LFC threshold for 'Down-regulated' status.
#' @param sig Numeric. The significance threshold (e.g., 0.05) to be used with
#'   `statval`.
#' @param adjval String. The p-value adjustment method for `p.adjust`
#'   (e.g., "BH").
#' @param statval Numeric. Which p-value to use for significance: `1` for raw
#'   p-value, `2` for adjusted p-value.
#' @param unique_proteins A list with two elements containing proteins unique to
#'   the control (`[1]`) and treatment (`[2]`) groups. Used only if `way = 2`.
#' @param way Numeric. `1` to analyze only shared proteins, `2` to include
#'   unique proteins from `unique_proteins`.
#' @param psms Logical. If `test = 2` (limma), `TRUE` to use `DEqMS` for
#'   spectral count weighting.
#' @param platform Numeric. Used if `psms = TRUE`. Specifies data source for
#'   finding PSM columns: `1`=MaxQuant, `2`=MSFragger, `4`=Proteome Discoverer.
#' @param selected_conditions A character vector of length 2.
#'   `selected_conditions[2]` is treated as the control group, and
#'   `selected_conditions[1]` as the treatment group.
#' @param diann_dir String. Directory path for DIA-NN. **Note:** This parameter
#'   is currently not implemented in the function body.
#'
#' @return A data frame (`results.eb`) with DE results, including `logFC`,
#'   `p.value`, `adj.P.Val`, `expression`, `Protein`, and
#'   `Protein_description`. Additional columns may be present for `limma`/`DEqMS`
#'   (e.g., `sca.P.Value`).
#'
#' @details
#' This function acts as a wrapper for various statistical tests. The `logFC` is
#' calculated as `mean(treatment) - mean(control)`.
#'
#' When `way = 2`, unique proteins are added to the results with imputed `logFC`
#' and p-values (e.g., p-value = min(p.value)/100) to ensure they appear in
#' visualizations.
#'
#' The function uses several non-standard variable assignment methods
#' (e.g., `assign`, `ls`, `get`) to dynamically create sample groups.
#'
#' @importFrom dplyr filter pull case_when
#' @importFrom limma lmFit eBayes topTable
#' @importFrom DEqMS spectraCounteBayes
#' @importFrom stats model.matrix p.adjust t.test wilcox.test
#' @importFrom matrixStats rowMins
#'
#' @export
#'
#' @examples
#' # 1. Create a mock data frame (df)
#' set.seed(123)
#' df_data <- data.frame(
#'   Protein = paste0("Prot", 1:5),
#'   Protein_description = paste("Description for", paste0("Prot", 1:5)),
#'   Cond1_R1 = rnorm(5, 10, 1),
#'   Cond1_R2 = rnorm(5, 10, 1),
#'   Cond2_R1 = rnorm(5, 12, 1), # Cond2 is higher
#'   Cond2_R2 = rnorm(5, 12, 1)
#' )
#' rownames(df_data) <- df_data$Protein
#'
#' # 2. Create mock metadata
#' meta <- data.frame(
#'   log2_col = c("Cond1_R1", "Cond1_R2", "Cond2_R1", "Cond2_R2"),
#'   group = c("Condition1", "Condition1", "Condition2", "Condition2")
#' )
#'
#' # 3. Define unique proteins (for way = 2)
#' uniques_list <- list(
#'   data.frame(Protein = "Prot_Control", Protein_description = "Only in Control"),
#'   data.frame(Protein = "Prot_Treat", Protein_description = "Only in Treatment")
#' )
#'
#' # Run analysis (t-test, unpaired, no PSMs, include uniques)
#' \dontrun{
#' results <- statistical_analysis(
#'   df = df_data,
#'   test = 1, # t-test
#'   paired = FALSE,
#'   metadata = meta,
#'   logfcup = 1,
#'   logfcdown = -1,
#'   sig = 0.05,
#'   adjval = "BH",
#'   statval = 2, # Use adjusted p-value
#'   unique_proteins = uniques_list,
#'   way = 2, # Include uniques
#'   psms = FALSE,
#'   platform = 1,
#'   selected_conditions = c("Condition2", "Condition1") # Treat 2 vs 1
#' )
#'
#' print(results)
#' }
#'
statistical_analysis <- function(df, test, paired = FALSE, metadata, logfcup, logfcdown, sig, adjval, statval, unique_proteins, way, psms, platform, selected_conditions, diann_dir = NULL){

  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]

  condition1_names <- metadata %>%
    dplyr::filter(group == first_group) %>%
    dplyr::pull(log2_col)

  replicas_condicion1 <- length(condition1_names)

  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]

  condition2_names <- metadata %>%
    dplyr::filter(group == second_group) %>%
    dplyr::pull(log2_col)

  replicas_condicion2 <- length(condition2_names)

  condi.names <- c(condition1_names,condition2_names)

  if (test == 2){
    if (psms == TRUE){
      #Control columns
      for (i in 1:replicas_condicion1){
        nam <- paste("control_sample", i, sep = "")
        assign(nam, condition1_names[i])
      }

      #Treatment columns
      for (i in 1:replicas_condicion2){
        nam <- paste("prob_column", i, sep = "")
        assign(nam, condition2_names[i])
      }
      ct <- c()
      for (i in ls()[grep("control_sample", ls())]){
        new_value_control <- get(i)
        ct <- c(ct, new_value_control)
      }
      print(ct)

      tr <- c()
      for (i in ls()[grep("prob_column", ls())]){
        new_value_prob <- get(i)
        tr <- c(tr, new_value_prob)
      }
      print(tr)

      if (paired == FALSE){
        control <- rep(1, replicas_condicion1)
        treatment <- rep(2, replicas_condicion2)
        design <- model.matrix(~factor(c(control, treatment))) # base R (stats)
      } else if (paired == TRUE){
        pairinfo = factor(rep(1:replicas_condicion1,2))
        control <- rep(1, replicas_condicion1)
        treatment <- rep(2, replicas_condicion2)
        design <- model.matrix(~pairinfo+factor(c(control, treatment))) # base R (stats)
      }

      dat <- df[, c(ct, tr)]
      n <- dim(dat)[1]
      fit <- limma::lmFit(dat, design)
      fit.eb <- limma::eBayes(fit)
      if (platform == 1){ #MaxQuant
        count_columns = grep("Razor...unique.peptides", colnames(df))
        df[count_columns] <- sapply(df[count_columns], as.numeric)
        psm.count.table = data.frame(count = matrixStats::rowMins(
          as.matrix(df[,count_columns])), row.names =  df$Protein)

        #psm.count.table[psm.count.table$count == 0] <- 1 # Añadimos solo un uno para los casos en los que el valor de PSMs es igual a 0
        psm.count.table$count[psm.count.table$count == 0] <- 1
        fit.eb$count <- psm.count.table$count

      } else if (platform == 2){ #MSFragger
        count_columns <- grep("Combined.Unique.Spectral.Count", colnames(df))
        count_column <- as.numeric(df$Combined.Unique.Spectral.Count)
        count_column <- replace(count_column, count_column==0, 1)
        fit.eb$count <- count_column

      } else if (platform == 4){ #Proteome Discoverer
        count_columns <- as.numeric(df$`# PSMs`)
        fit.eb$count <- count_columns

      }
      fit <- DEqMS::spectraCounteBayes(fit.eb)

      coef_col <- 2 #Para acceder a la columna de LogFC y P value

      results.eb = limma::topTable(fit,coef = coef_col,n= Inf)
      p.value.column <- grep("P.Value", colnames(results.eb))
      colnames(results.eb)[p.value.column] <- "p.value"
      results.eb$adj.P.Val  = p.adjust(results.eb$p.value, # base R (stats)
                                       method = adjval)

      results.eb$Protein = as.numeric(rownames(results.eb))
      #results.eb$count = fit$count[results.eb$Protein]

      results.eb$sca.t = fit$sca.t[results.eb$Protein,coef_col]
      results.eb$sca.P.Value = as.numeric(fit$sca.p[results.eb$Protein,coef_col])
      results.eb$sca.adj.pval = as.numeric(p.adjust(results.eb$sca.P.Value, # base R (stats)
                                                    method = "BH"))
      results.eb = results.eb[order(results.eb$sca.P.Value), ]

      if (statval == 1){
        expression <- dplyr::case_when(results.eb$logFC >= logfcup & -log10(results.eb$sca.P.Value) >= -log10(sig) ~ "Up-regulated",
                                       results.eb$logFC <= logfcdown & -log10(results.eb$sca.P.Value) >= -log10(sig) ~ "Down-regulated",
                                       TRUE ~ "Unchanged")#labels para expresion
      } else if (statval == 2){
        expression <- dplyr::case_when(results.eb$logFC >= logfcup & -log10(results.eb$sca.adj.pval) >= -log10(sig) ~ "Up-regulated",
                                       results.eb$logFC <= logfcdown & -log10(results.eb$sca.adj.pval) >= -log10(sig) ~ "Down-regulated",
                                       TRUE ~ "Unchanged")
      }
      results.eb$expression <- expression
      results_rownames <- rownames(results.eb)
      results.eb$Protein <- df[results_rownames, "Protein"]
      results.eb$Protein_description <- df[results_rownames, "Protein_description"]

      row.names(results.eb) <- results.eb$Protein
      #row.names(results.eb) <- make.unique(results.eb$Protein)
      results.eb <- results.eb[-c(2,3,6,8)]
      results.eb <- results.eb[c(1,2,3,4,8,6,5,7)]
      if (way == 1){
        return(results.eb)
      } else if (way == 2){
        #Empezamos a trabajar con las unicas control
        unique.control <- as.data.frame(unique_proteins[1])
        fusioncontrol.df <- data.frame(matrix(ncol = 8))
        n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description", "sca.P.Value", "sca.adj.pval")
        colnames(fusioncontrol.df) <- n

        for (i in unique.control$Protein){
          new <- c(min(results.eb$logFC, na.rm = TRUE)-2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Down-regulated", i, unique.control$Protein_description[unique.control$Protein==i], min(results.eb$sca.P.Value, na.rm = TRUE)/100, min(results.eb$sca.adj.pval, na.rm = TRUE)/100)
          fusioncontrol.df <- rbind(fusioncontrol.df, new)
        }
        fusioncontrol.df <- fusioncontrol.df[-1,]

        unique.treatment <- as.data.frame(unique_proteins[2])
        fusiontreatment.df <- data.frame(matrix(ncol = 8))
        n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description", "sca.P.Value", "sca.adj.pval")
        colnames(fusiontreatment.df) <- n

        for (i in unique.treatment$Protein){
          new <- c(max(results.eb$logFC, na.rm = TRUE)+2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Up-regulated", i, unique.treatment$Protein_description[unique.treatment$Protein==i], min(results.eb$sca.P.Value, na.rm = TRUE)/100, min(results.eb$sca.adj.pval, na.rm = TRUE)/100)
          fusiontreatment.df <- rbind(fusiontreatment.df, new)
        }
        fusiontreatment.df <- fusiontreatment.df[-1,]


        unique.proteins.limma <- rbind(fusioncontrol.df, fusiontreatment.df)
        row.names(unique.proteins.limma) <- unique.proteins.limma$Protein


        unique.proteins.limma[c("logFC", "p.value", "adj.P.Val", "sca.P.Value", "sca.adj.pval")] <- sapply(unique.proteins.limma[c("logFC", "p.value", "adj.P.Val", "sca.P.Value", "sca.adj.pval")], as.numeric)

        results.eb <- rbind(unique.proteins.limma, results.eb)
        #row.names(results.eb) <- make.unique(results.eb$Protein)
        return(results.eb)
      }

    } else if (psms == FALSE){
      #Control columns
      for (i in 1:replicas_condicion1){
        nam <- paste("control_sample", i, sep = "")
        assign(nam, condition1_names[i])
      }

      #Treatment columns
      for (i in 1:replicas_condicion2){
        nam <- paste("prob_column", i, sep = "")
        assign(nam, condition2_names[i])
      }
      ct <- c()
      for (i in ls()[grep("control_sample", ls())]){
        new_value_control <- get(i)
        ct <- c(ct, new_value_control)
      }
      print(ct)

      tr <- c()
      for (i in ls()[grep("prob_column", ls())]){
        new_value_prob <- get(i)
        tr <- c(tr, new_value_prob)
      }
      print(tr)

      if (paired == FALSE){
        control <- rep(1, replicas_condicion1)
        treatment <- rep(2, replicas_condicion2)
        design <- model.matrix(~factor(c(control, treatment))) # base R (stats)
      } else if (paired == TRUE){
        pairinfo = factor(rep(1:replicas_condicion1,2))
        control <- rep(1, replicas_condicion1)
        treatment <- rep(2, replicas_condicion2)
        design <- model.matrix(~pairinfo+factor(c(control, treatment))) # base R (stats)
      }

      dat <- df[, c(ct, tr)]
      n <- dim(dat)[1]
      fit <- limma::lmFit(dat, design)
      fit.eb <- limma::eBayes(fit)
      print(colnames(fit.eb))
      logFC <- fit.eb$coefficients[, 2] #Calculo del log fold-change
      p.value <- fit.eb$p.value[, 2]    # p-valor moderado correspondiente al estad?stico t moderado.
      adj.P.Val <- p.adjust(p.value, method = adjval) # base R (stats)

      if (statval == 1){
        expression <- dplyr::case_when(logFC >= logfcup & -log10(p.value) >= -log10(sig) ~ "Up-regulated",
                                       logFC <= logfcdown & -log10(p.value) >= -log10(sig) ~ "Down-regulated",
                                       TRUE ~ "Unchanged")#labels para expresion
      } else if (statval == 2){
        expression <- dplyr::case_when(logFC >= logfcup & -log10(adj.P.Val) >= -log10(sig) ~ "Up-regulated",
                                       logFC <= logfcdown & -log10(adj.P.Val) >= -log10(sig) ~ "Down-regulated",
                                       TRUE ~ "Unchanged")
      }
      results.eb <- data.frame(logFC, p.value, adj.P.Val, expression)
      results_rownames <- rownames(results.eb) #Obtenemos el nombre de las filas

      # Extract Protein and Protein_description directly using vectorized operations
      results.eb$Protein <- df[results_rownames, "Protein"]
      results.eb$Protein_description <- df[results_rownames, "Protein_description"]
      #row.names(results.eb) <- make.unique(results.eb$Protein)

      if (way == 1){
        return(results.eb)
      } else if (way == 2){
        #Empezamos a trabajar con las unicas control
        unique.control <- as.data.frame(unique_proteins[1])
        fusioncontrol.df <- data.frame(matrix(ncol = 6))
        n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
        colnames(fusioncontrol.df) <- n

        for (i in unique.control$Protein){
          new <- c(min(results.eb$logFC, na.rm = TRUE)-2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Down-regulated", i, unique.control$Protein_description[unique.control$Protein==i] )
          fusioncontrol.df <- rbind(fusioncontrol.df, new)
        }
        fusioncontrol.df <- fusioncontrol.df[-1,]

        #Empezamos a trabajar con la unicas tratamiento
        unique.treatment <- as.data.frame(unique_proteins[2])
        fusiontreatment.df <- data.frame(matrix(ncol = 6))
        n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
        colnames(fusiontreatment.df) <- n

        for (i in unique.treatment$Protein){
          new <- c(max(results.eb$logFC, na.rm = TRUE)+2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Up-regulated", i, unique.treatment$Protein_description[unique.treatment$Protein==i] )
          fusiontreatment.df <- rbind(fusiontreatment.df, new)
        }
        fusiontreatment.df <- fusiontreatment.df[-1,]


        unique.proteins.limma <- rbind(fusioncontrol.df, fusiontreatment.df)
        row.names(unique.proteins.limma) <- unique.proteins.limma$Protein


        unique.proteins.limma[c("logFC", "p.value", "adj.P.Val")] <- sapply(unique.proteins.limma[c("logFC", "p.value", "adj.P.Val")], as.numeric)

        results.eb <- rbind(unique.proteins.limma, results.eb)
        #row.names(results.eb) <- make.unique(results.eb$Protein)
        return(results.eb)
      }
    }

  } else if (test == 1){

    valuesA_cols <- df[, condition1_names]
    valuesB_cols <- df[, condition2_names]

    for (i in unique(seq_len(nrow(df)))) {
      valuesA <- as.numeric(valuesA_cols[i, ])
      valuesB <- as.numeric(valuesB_cols[i, ])

      if (sum(is.finite(valuesA)) >= 2 && sum(is.finite(valuesB)) >= 2) {
        testResults <- t.test(x = valuesA, y = valuesB, paired = paired) # base R (stats)

        df$pValue[i] <- testResults$p.value
        df$tStat[i] <- testResults$statistic
        df$logFC[i] <- mean(valuesB, na.rm = TRUE) - mean(valuesA, na.rm = TRUE)

        # Storing results in a separate data frame
        results.eb <- data.frame(df$logFC, df$pValue)
      } else {
        df$pValue[i] <- NA
        df$tStat[i] <- NA
        df$logFC[i] <- NA
      }
    }


  } else if (test == 3){
    valuesA_cols <- df[, condition1_names]
    valuesB_cols <- df[, condition2_names]

    for (i in unique(seq_len(nrow(df)))) {
      valuesA <- as.numeric(valuesA_cols[i, ])
      valuesB <- as.numeric(valuesB_cols[i, ])

      if (sum(is.finite(valuesA)) >= 2 && sum(is.finite(valuesB)) >= 2) {
        testResults <- wilcox.test(x = valuesA, y = valuesB, paired = paired) # base R (stats)

        df$pValue[i] <- testResults$p.value
        df$tStat[i] <- testResults$statistic
        df$logFC[i] <- mean(valuesB, na.rm = TRUE) - mean(valuesA, na.rm = TRUE)

        # Storing results in a separate data frame
        results.eb <- data.frame(df$logFC, df$pValue)
      } else {
        df$pValue[i] <- NA
        df$tStat[i] <- NA
        df$logFC[i] <- NA
      }
    }


  }
  colnames(results.eb) <- c("logFC", "p.value")
  results.eb$adj.P.Val <- p.adjust(results.eb$p.value, method = adjval) # base R (stats)
  if (statval==1){
    results.eb$expression <- dplyr::case_when(results.eb$logFC >= logfcup & -log10(results.eb$p.value) >= -log10(sig) ~ "Up-regulated",
                                              results.eb$logFC <= logfcdown & -log10(results.eb$p.value) >= -log10(sig) ~ "Down-regulated",
                                              TRUE ~ "Unchanged")#labels para expresion
  } else if (statval==2){
    results.eb$expression <- dplyr::case_when(results.eb$logFC >= logfcup & -log10(results.eb$adj.P.Val) >= -log10(sig) ~ "Up-regulated",
                                              results.eb$logFC <= logfcdown & -log10(results.eb$adj.P.Val) >= -log10(sig) ~ "Down-regulated",
                                              TRUE ~ "Unchanged")
  }


  # Assuming results.eb already exists and has the appropriate structure
  results_rownames <- rownames(results.eb)

  # Extract Protein and Protein_description directly using vectorized operations
  results.eb$Protein <- df[results_rownames, "Protein"]
  results.eb$Protein_description <- df[results_rownames, "Protein_description"]

  #row.names(results.eb) <- make.unique(results.eb$Protein)

  if (way == 1){
    return(results.eb)
  } else if (way == 2){
    #Empezamos a trabajar con las unicas control
    unique.control <- as.data.frame(unique_proteins[1])
    fusioncontrol.df <- data.frame(matrix(ncol = 6))
    n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
    colnames(fusioncontrol.df) <- n

    for (i in unique.control$Protein){
      new <- c(min(results.eb$logFC, na.rm = TRUE)-2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Down-regulated", i, unique.control$Protein_description[unique.control$Protein==i] )
      fusioncontrol.df <- rbind(fusioncontrol.df, new)
    }

    fusioncontrol.df <- fusioncontrol.df[-1,]

    #Empezamos a trabajar con la unicas tratamiento
    unique.treatment <- as.data.frame(unique_proteins[2])
    fusiontreatment.df <- data.frame(matrix(ncol = 6))
    n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
    colnames(fusiontreatment.df) <- n

    for (i in unique.treatment$Protein){
      new <- c(max(results.eb$logFC, na.rm = TRUE)+2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Up-regulated", i, unique.treatment$Protein_description[unique.treatment$Protein==i] )
      fusiontreatment.df <- rbind(fusiontreatment.df, new)
    }
    fusiontreatment.df <- fusiontreatment.df[-1,]


    unique.proteins.limm <- rbind(fusioncontrol.df, fusiontreatment.df)
    row.names(unique.proteins.limm) <- unique.proteins.limm$Protein

    unique.proteins.limm[c("logFC", "p.value", "adj.P.Val")] <- sapply(unique.proteins.limm[c("logFC", "p.value", "adj.P.Val")], as.numeric)

    results.eb <- rbind(unique.proteins.limm, results.eb)
    #row.names(results.eb) <- make.unique(results.eb$Protein)
    return(results.eb)
  }
}


##################################################################################
#Plots
#Create an Interactive Volcano Plot
#'
#' This function generates an interactive volcano plot using the 'plotly' package
#' based on differential expression results.
#'
#' @param limma A data frame containing differential expression results.
#'   This data frame must include columns:
#'   \itemize{
#'     \item \code{logFC}: Log-fold change values.
#'     \item \code{Protein}: Protein identifiers (used for hover-text).
#'     \item \code{expression}: A character vector (e.g., "Up-regulated",
#'       "Down-regulated", "Unchanged") used for coloring.
#'     \item \code{p.value} or \code{sca.P.Value}: Nominal p-values.
#'     \item \code{adj.P.Val} or \code{sca.adj.pval}: Adjusted p-values.
#'   }
#' @param title A character string for the plot title.
#' @param label A numeric value specifying the size of the markers (points)
#'   in the plot. Passed to `size = I(label)`.
#' @param statval A numeric value (1 or 2) specifying which p-value to plot on
#'   the y-axis.
#'   \itemize{
#'     \item \code{1}: Plots the nominal p-value (e.g., `p.value` or
#'       `sca.P.Value`).
#'     \item \code{2}: Plots the adjusted p-value (e.g., `adj.P.Val` or
#'       `sca.adj.pval`).
#'   }
#' @param psms A logical value (`TRUE` or `FALSE`).
#'   \itemize{
#'     \item \code{TRUE}: Uses PSM-aware p-values from DEqMS (i.e.,
#'       `sca.P.Value` or `sca.adj.pval`).
#'     \item \code{FALSE}: Uses standard p-values (i.e., `p.value` or
#'       `adj.P.Val`).
#'   }
#'
#' @return A `plotly` object, which renders as an interactive volcano plot.
#'
#' @details
#' The function selects the appropriate p-value column based on the `statval`
#' and `psms` arguments. The y-axis label is automatically set to "-log10 p-value"
#' or "-log10 q-value" (for adjusted p-values). Points are colored based on
#' the `expression` column.
#'
#' @examples
#' \dontrun{
#' # Create a mock data frame for demonstration
#' limma_results <- data.frame(
#'   logFC = rnorm(100, 0, 2),
#'   Protein = paste0("PROT", 1:100),
#'   p.value = runif(100, 0, 1),
#'   adj.P.Val = runif(100, 0, 1),
#'   sca.P.Value = runif(100, 0, 1),
#'   sca.adj.pval = runif(100, 0, 1),
#'   expression = sample(c("Up-regulated", "Down-regulated", "Unchanged"),
#'                       100, replace = TRUE)
#' )
#'
#' # Generate an interactive plot using adjusted p-values without PSM weighting
#' volcano_plot(limma = limma_results,
#'              title = "Interactive Volcano Plot",
#'              label = 5,
#'              statval = 2,
#'              psms = FALSE)
#' }
#'
#' @importFrom plotly plot_ly layout
#' @export

volcano_plot <- function(limma, title, label, statval, psms){

  logFC <- limma$logFC
  protein_ids <- limma$Protein

  #top <- readline("Introduzca cuantas proteinas quiere etiquetar:")
  if (statval == 1){

    if (psms == TRUE){
      plot <- plotly::plot_ly(data = limma, x = ~logFC, y = ~-log10(sca.P.Value), text = protein_ids,
                              type = "scatter", mode = "markers",
                              color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                              size = I(label))
    } else if( psms == FALSE){
      plot <- plotly::plot_ly(data = limma, x = ~logFC, y = ~-log10(p.value), text = protein_ids,
                              type = "scatter", mode = "markers",
                              color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                              size = I(label))
    }
    plot <- plot %>% plotly::layout(title = title,
                                    xaxis = list(title = list(text ='Log2 Fold Change')),
                                    yaxis = list(title = list(text = 'Log10 p-value')))

    return(plot)
  } else if (statval == 2){

    logFC <- limma$logFC
    protein_ids <- limma$Protein

    if (psms == TRUE){
      plot <- plotly::plot_ly(data = limma, x = ~logFC, y = ~-log10(sca.adj.pval), text = protein_ids,
                              type = "scatter", mode = "markers",
                              color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                              size = I(label))
    } else if( psms == FALSE){
      plot <- plotly::plot_ly(data = limma, x = ~logFC, y = ~-log10(adj.P.Val), text = protein_ids,
                              type = "scatter", mode = "markers",
                              color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                              size = I(label))
    }
    plot <- plot %>% plotly::layout(title = title,
                                    xaxis = list(title = list(text ='Log2 Fold Change')),
                                    yaxis = list(title = list(text = 'Log10 q-value')))

    return(plot)
  }

}

#' Create a Static Volcano Plot for Publication
#'
#' This function generates a static volcano plot using the 'ggplot2' package,
#' suitable for saving as a high-resolution image (e.g., TIFF) for
#' publication.
#'
#' @param limma A data frame containing differential expression results.
#'   This data frame must include columns:
#'   \itemize{
#'     \item \code{logFC} or \code{log2FC}: Log-fold change values.
#'     \item \code{Protein}: Protein identifiers.
#'     \item \code{expression}: A character vector (e.g., "Up-regulated",
#'       "Down-regulated", "Unchanged") used for coloring.
#'     \item \code{p.value} or \code{sca.P.Value}: Nominal p-values.
#'     \item \code{adj.P.Val} or \code{sca.adj.pval}: Adjusted p-values.
#'   }
#' @param title A character string for the plot title.
#' @param label A numeric value specifying the size of the points.
#'   \strong{Note:} In the current function, this parameter is only used when
#'   `statval = 1` and `psms = TRUE`. In all other cases, the size is
#'   hard-coded to `2`.
#' @param statval A numeric value (1 or 2) specifying which p-value to plot on
#'   the y-axis.
#'   \itemize{
#'     \item \code{1}: Plots the nominal p-value (e.g., `p.value` or
#'       `sca.P.Value`).
#'     \item \code{2}: Plots the adjusted p-value (e.g., `adj.P.Val` or
#'       `sca.adj.pval`).
#'   }
#' @param psms A logical value (`TRUE` or `FALSE`).
#'   \itemize{
#'     \item \code{TRUE}: Uses PSM-aware p-values from DEqMS (i.e.,
#'       `sca.P.Value` or `sca.adj.pval`).
#'     \item \code{FALSE}: Uses standard p-values (i.e., `p.value` or
#'       `adj.P.Val`).
#'   }
#'
#' @return A `ggplot` object, which renders as a static volcano plot.
#'
#' @details
#' The function uses `ggplot2::theme_light()` and customizes axis text and
#' line sizes. The y-axis label is automatically set to "-log10 p-value"
#' or "-log10 q-value" (for adjusted p-values).
#'
#' @examples
#' \dontrun{
#' # Create a mock data frame for demonstration
#' limma_results <- data.frame(
#'   logFC = rnorm(100, 0, 2),
#'   log2FC = limma_results$logFC, # For compatibility
#'   Protein = paste0("PROT", 1:100),
#'   p.value = runif(100, 0, 1),
#'   adj.P.Val = runif(100, 0, 1),
#'   sca.P.Value = runif(100, 0, 1),
#'   sca.adj.pval = runif(100, 0, 1),
#'   expression = sample(c("Up-regulated", "Down-regulated", "Unchanged"),
#'                       100, replace = TRUE)
#' )
#'
#' # Generate a static plot
#' static_plot <- volcano_plot_tiff(limma = limma_results,
#'                                  title = "Static Volcano Plot",
#'                                  label = 3,
#'                                  statval = 2,
#'                                  psms = TRUE)
#'
#' # You can then save it, e.g., with ggsave()
#' # ggplot2::ggsave("my_volcano_plot.tiff", static_plot,
#' #                 width = 7, height = 5, dpi = 300)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual guides
#'   guide_legend ggtitle labs theme theme_light element_text element_line
#' @export
volcano_plot_tiff <- function(limma, title, label,  statval, psms){

  plot <- NULL  # Initialize

  if (statval == 1){

    if(psms == TRUE){
      logFC <- limma$log2FC
      protein_ids <- limma$Protein

      plot <- ggplot2::ggplot(limma, ggplot2::aes(logFC, -log10(sca.P.Value))) +
        ggplot2::theme_light() +
        ggplot2::geom_point(ggplot2::aes(color = expression), size = label) +
        ggplot2::scale_color_manual(values = c("green3", "gray78", "firebrick3")) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
        ggplot2::ggtitle(title) +
        ggplot2::labs(y = "-log10 p-value", x = "log2 Fold change") +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 14),
          line = ggplot2::element_line(linewidth = 1.75)
        )

    } else if(psms == FALSE){
      plot <- ggplot2::ggplot(limma, ggplot2::aes(logFC, -log10(p.value))) +
        ggplot2::theme_light() +
        ggplot2::geom_point(ggplot2::aes(color = expression), size = 2) +
        ggplot2::scale_color_manual(values = c("green3", "gray78", "firebrick3")) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
        ggplot2::ggtitle(title) +
        ggplot2::labs(y = "-log10 p-value", x = "log2 Fold change") +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 14),
          line = ggplot2::element_line(linewidth = 1.75)
        )
    }

    return(plot)
  } else if (statval == 2){

    logFC <- limma$log2FC
    protein_ids <- limma$Protein

    if(psms == TRUE){
      plot <- ggplot2::ggplot(limma, ggplot2::aes(logFC, -log10(sca.adj.pval))) +
        ggplot2::theme_light() +
        ggplot2::geom_point(ggplot2::aes(color = expression), size = 2) +
        ggplot2::scale_color_manual(values = c("green3", "gray78", "firebrick3")) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
        ggplot2::ggtitle(title) +
        ggplot2::labs(y = "-log10 q-value", x = "log2 Fold change") +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 14),
          line = ggplot2::element_line(linewidth = 1.75)
        )

    } else if(psms == FALSE){
      plot <- ggplot2::ggplot(limma, ggplot2::aes(logFC, -log10(adj.P.Val))) +
        ggplot2::theme_light() +
        ggplot2::geom_point(ggplot2::aes(color = expression), size = 2) +
        ggplot2::scale_color_manual(values = c("green3", "gray78", "firebrick3")) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.5))) +
        ggplot2::ggtitle(title) +
        ggplot2::labs(y = "-log10 q-value", x = "log2 Fold change") +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 14),
          line = ggplot2::element_line(linewidth = 1.75)
        )
    }
  }

  return(plot)
}

#' @title Custom PCA Plot (Base R)
#' @description Performs PCA and generates a 2D scatter plot of the first two
#' principal components using base R plotting.
#'
#' @param x A data frame or matrix with expression data. Columns are samples
#'          (e.g., LOG2 intensity values). Rows are features (e.g. proteins).
#' @param metadata A data frame with sample metadata. Must have 'group' and
#'                 'log2_col' columns.
#' @param selected_conditions A character vector of length 2 specifying the two
#'                            groups from the 'group' column to compare.
#'
#' @return This function is called for its side effect of generating a base R
#'         plot. It does not return an object.
#' @export
#' @importFrom dplyr %>% filter pull
#' @importFrom stats prcomp na.omit
#' @importFrom graphics plot text
#'
#' @examples
#' # Create synthetic data
#' x_data <- data.frame(
#'   LOG2.C1 = rnorm(10, 8), LOG2.C2 = rnorm(10, 8),
#'   LOG2.T1 = rnorm(10, 9), LOG2.T2 = rnorm(10, 9),
#'   row.names = paste0("Prot", 1:10)
#' )
#' metadata <- data.frame(
#'   log2_col = c("LOG2.C1", "LOG2.C2", "LOG2.T1", "LOG2.T2"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' selected_cond <- c("Treatment", "Control") # Matches user's logic [2] then [1]
#'
#' # Generate the plot
#' pca(x_data, metadata, selected_cond)
#'

pca <- function(x, metadata, selected_conditions, pc_x = 1, pc_y = 2) {

  # Convert character inputs to numeric for indexing
  pc_x <- as.numeric(pc_x)
  pc_y <- as.numeric(pc_y)

  # --- Sample Filtering (No change) ---
  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]
  condition1_names <- metadata %>%
    dplyr::filter(group == first_group) %>%
    dplyr::pull(log2_col)

  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]
  condition2_names <- metadata %>%
    dplyr::filter(group == second_group) %>%
    dplyr::pull(log2_col)

  group_names <- c(condition1_names, condition2_names)

  x_filtered <- na.omit(x[group_names]) # na.omit is base R
  tdata <- t(x_filtered) # t is base R
  pc <- prcomp(tdata) # prcomp is base R

  # --- Eigenvalue/Contribution Calculation (Modified for selected PCs) ---
  eigenValues <- pc$sdev^2
  totalEigenValue <- sum(eigenValues)

  # Check if selected PCs exist
  if (pc_x > length(eigenValues) || pc_y > length(eigenValues)) {
    stop("Selected principal component index is out of bounds.")
  }

  # Calculate contribution for the *selected* PCs
  contribution_x <- round(100 * eigenValues[pc_x] / totalEigenValue)
  contribution_y <- round(100 * eigenValues[pc_y] / totalEigenValue)

  # --- Data for Plotting (Modified to use selected PCs) ---
  pca_data <- data.frame(
    PC_X = pc$x[, pc_x], # Use pc_x for the x-axis
    PC_Y = pc$x[, pc_y], # Use pc_y for the y-axis
    Sample = rownames(pc$x),
    Condition = rep(c(first_group, second_group), times = c(length(condition1_names), length(condition2_names)))
  )

  # --- Plotting (Modified to use generic labels and selected contributions) ---
  plot_object <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC_X, y = PC_Y, color = Condition, label = Sample)) +
    ggplot2::geom_point(size = 4) +

    ggrepel::geom_text_repel() +
    ggplot2::labs(
      title = "Principal Component Analysis",
      # Use the selected PC numbers in the labels
      x = paste0("PC", pc_x, " [", contribution_x, "%]"),
      y = paste0("PC", pc_y, " [", contribution_y, "%]"),
      color = "Condition" # Legend title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(plot_object)
}
#' @title t-SNE Plot
#' @description Generates a t-SNE plot for sample clustering using ggplot2.
#'
#' @param x A data frame or matrix with expression data. Columns are samples
#'          (e.g., LOG2 intensity values). Rows are features (e.g. proteins).
#' @param metadata A data frame with sample metadata. Must have 'group' and
#'                 'log2_col' columns.
#' @param perplexity_num Numeric. The perplexity value for the t-SNE algorithm.
#' @param selected_conditions A character vector of length 2 specifying the two
#'                            groups from the 'group' column to compare.
#'
#' @return A ggplot object.
#' @export
#' @importFrom dplyr %>% filter pull
#' @importFrom stats na.omit
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot aes geom_point geom_text stat_ellipse labs
#' @importFrom ggplot2 theme_minimal
#'
#' @examples
#' # Create synthetic data
#' x_data <- data.frame(
#'   LOG2.C1 = rnorm(10, 8), LOG2.C2 = rnorm(10, 8),
#'   LOG2.T1 = rnorm(10, 9), LOG2.T2 = rnorm(10, 9),
#'   row.names = paste0("Prot", 1:10)
#' )
#' metadata <- data.frame(
#'   log2_col = c("LOG2.C1", "LOG2.C2", "LOG2.T1", "LOG2.T2"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' selected_cond <- c("Treatment", "Control")
#'
#' # Rtsne requires perplexity < (n-1)/3. Here n=4, so perplexity must be 1.
#' if (requireNamespace("Rtsne")) {
#'   tsne(x_data, metadata, 1, selected_cond)
#' }
#'
tsne <- function(x, metadata, perplexity_num, selected_conditions){

  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]

  condition1_names <- metadata %>%
    filter(group == first_group) %>%
    pull(log2_col)

  rep1 <- length(condition1_names)

  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]

  condition2_names <- metadata %>%
    filter(group == second_group) %>%
    pull(log2_col)

  rep2 <- length(condition2_names)

  group <- c(condition1_names,condition2_names)
  x <- na.omit(x[group])

  tdata<-as.data.frame(t(x[group]))

 #For reproducibility
  tsne <- Rtsne(X = tdata, perplexity = as.numeric(perplexity_num))

  sample_names <- rownames(tdata)

  my_colors <-  rep(c("salmon", "blue"), c(rep1, rep2))

  tsne_df <- data.frame(
    Dim1 = tsne$Y[, 1],
    Dim2 = tsne$Y[, 2],
    Sample = sample_names,
    Group = rep(c(first_group, second_group), c(rep1, rep2)) #Mejorar como definimos los grupos
  )

  ggplot(tsne_df, aes(x = Dim1, y = Dim2, label = Sample, color = Group)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    stat_ellipse(aes(group = Group), linetype = "dashed", size = 1) +
    labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal()

}

#' @title Base R Heatmap
#' @description Creates a heatmap using the base R `heatmap` function.
#'
#' @param data A data frame containing expression data and a 'Protein' column.
#' @param cond.names A character vector of the column names (samples) to include
#'                   in the heatmap.
#' @param title A character string for the heatmap's title.
#'
#' @return Returns (invisibly) a list with components from the `heatmap`
#'         function, such as `rowInd` and `colInd`.
#' @export
#' @importFrom stats heatmap
#' @importFrom gplots greenred
#'
#' @examples
#' # Create synthetic data
#' heat_data <- data.frame(
#'   Protein = paste0("Prot", 1:5),
#'   Sample1 = rnorm(5, 10),
#'   Sample2 = rnorm(5, 12)
#' )
#' sample_names <- c("Sample1", "Sample2")
#'
#' # Generate the plot
#' my_heatmap(heat_data, sample_names, "My Heatmap")
#'

my_heatmap <- function(data, cond.names, title){

  heatmap(as.matrix(data[cond.names]), labRow = data$Protein, main = title, Rowv = NULL, Colv = NA, col = gplots::greenred(75), cexCol = 0.6)



}
#' @title Differential Protein Heatmap (ggplot)
#' @description Creates a ggplot2 heatmap of differentially expressed proteins.
#'
#' @param limma A data frame from statistical analysis, must have 'Protein' and
#'              'expression' columns (e.g., 'Up-regulated', 'Down-regulated').
#' @param data A data frame containing expression data, must have 'Protein' and
#'             sample columns (e.g., LOG2 intensity values).
#' @param cond.names A character vector of the column names (samples) to include
#'                   in the heatmap.
#' @param title A character string for the heatmap's title.
#'
#' @return A ggplot object.
#' @export
#' @importFrom dplyr %>% bind_rows filter all_of group_by mutate ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal
#' @importFrom ggplot2 theme element_text annotate theme_void
#'
#' @examples
#' # Create synthetic data
#' limma_results <- data.frame(
#'   Protein = c("Prot1", "Prot2", "Prot3"),
#'   expression = c("Up-regulated", "Not Significant", "Down-regulated")
#' )
#' data_expr <- data.frame(
#'   Protein = c("Prot1", "Prot2", "Prot3"),
#'   Sample1 = c(12, 10, 8),
#'   Sample2 = c(13, 10, 7)
#' )
#' sample_names <- c("Sample1", "Sample2")
#'
#' # Generate the plot
#' my_heatmap_differential(limma_results, data_expr, sample_names, "DE Heatmap")
#'
my_heatmap_differential <- function(limma, data, cond.names, title) {

  top_proteins <- dplyr::bind_rows(
    limma %>% dplyr::filter(expression == 'Up-regulated'),
    limma %>% dplyr::filter(expression == 'Down-regulated')
  )

  if (nrow(top_proteins) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x=0, y=0, label="No differentially expressed proteins found.") +
             ggplot2::theme_void())
  }

  heatmap_data_long <- data %>%
    dplyr::filter(Protein %in% top_proteins$Protein) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cond.names),
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    dplyr::mutate(Intensity = as.numeric(Intensity))

  heatmap_data_long <- heatmap_data_long %>%
    dplyr::group_by(Protein) %>%
    dplyr::mutate(Intensity_scaled = scale(Intensity)) %>% # scale is base R
    dplyr::ungroup()

  heatmap_object <- ggplot2::ggplot(heatmap_data_long, ggplot2::aes(x = Sample, y = Protein, fill = Intensity_scaled)) +
    ggplot2::geom_tile(color = "white") +

    ggplot2::scale_fill_gradient2(
      low = "green",
      mid = "black",
      high = "red",
      midpoint = 0, # Assumes scaled data is centered around 0
      name = "Scaled Intensity"
    ) +
    ggplot2::labs(title = title, x = "Samples", y = "Proteins") +
    ggplot2::theme_minimal() + # A clean theme for heatmaps
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(heatmap_object)
}


#' @title Differential Protein Boxplot
#' @description Creates a boxplot with jittered points for a single protein,
#' comparing two conditions.
#'
#' @param df A data frame with expression data. Must have a 'Protein' column
#'           and LOG2 intensity columns.
#' @param metadata A data frame with sample metadata. Must have 'group' and
#'                 'log2_col' columns.
#' @param protein A character string specifying the single protein (from the
#'                'Protein' column) to plot.
#' @param LOG2.names A character vector of all LOG2 sample column names to use.
#' @param selected_conditions A character vector of length 2 specifying the two
#'                            groups from the 'group' column to compare.
#'
#' @return This function is called for its side effect of printing a ggplot.
#'         It invisibly returns the ggplot object `p`.
#' @export
#' @importFrom dplyr %>% filter pull
#' @importFrom ggplot2 ggplot aes geom_jitter position_dodge geom_boxplot labs
#' @importFrom ggplot2 scale_color_brewer facet_wrap theme element_blank
#' @importFrom DOSE theme_dose
#'
#' @examples
#' # Create synthetic data
#' df_box <- data.frame(
#'   Protein = "Prot1",
#'   LOG2.C1 = rnorm(1, 8), LOG2.C2 = rnorm(1, 8),
#'   LOG2.T1 = rnorm(1, 9), LOG2.T2 = rnorm(1, 9)
#' )
#' metadata_box <- data.frame(
#'   log2_col = c("LOG2.C1", "LOG2.C2", "LOG2.T1", "LOG2.T2"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' log_names <- c("LOG2.C1", "LOG2.C2", "LOG2.T1", "LOG2.T2")
#' selected_cond <- c("Treatment", "Control")
#'
#' # Generate the plot
#' if (requireNamespace("DOSE")) {
#'   Diferential_boxplot(df_box, metadata_box, "Prot1", log_names, selected_cond)
#' }
#'

Diferential_boxplot <- function(df, metadata, protein, LOG2.names, selected_conditions){
  row.names(df) <- df$Protein

  subset_data <- df[df$Protein == protein, c(LOG2.names)]

  if (nrow(subset_data) == 0) {
    stop(paste("Protein", protein, "not found in df$Protein"))
  }

  first_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[2]]
  second_group <- unique(metadata$group)[unique(metadata$group) == selected_conditions[1]]

  cond.colums <- intersect(colnames(subset_data), metadata %>%
                             dplyr::filter(group == first_group) %>%
                             dplyr::pull(log2_col))

  trat.colums <- intersect(colnames(subset_data), metadata %>%
                             dplyr::filter(group == second_group) %>%
                             dplyr::pull(log2_col))

  cond.colums <- as.character(cond.colums)
  trat.colums <- as.character(trat.colums)

  all_samples <- c(cond.colums, trat.colums)

  dt <- data.frame(
    Condition = rep(c(first_group, second_group),
                    c(length(cond.colums), length(trat.colums))),
    val = as.vector(unlist(subset_data[, all_samples])),
    Samples = c(cond.colums, trat.colums)
  )

  p <- ggplot2::ggplot(dt, ggplot2::aes(Condition, val)) +
    DOSE::theme_dose() +
    ggplot2::geom_jitter(ggplot2::aes(color = factor(Samples)),
                         size = 5, position = ggplot2::position_dodge(width=0.4)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      y = expression(log[2]~"Intensity"),
      col = "Samples") +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::facet_wrap(~{protein}) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  print(p)
}



##################################################################################
#Functional analysis
#' @title Find GO Terms (gprofiler2)
#' @description Performs Gene Ontology (GO) and pathway enrichment analysis
#' for up- and down-regulated proteins using gprofiler2.
#'
#' @param df A data frame from statistical analysis (e.g., limma), must have
#'           'expression', 'p.value', 'Protein', and/or 'Protein_description'
#'           columns.
#' @param raw A data frame of the "raw" data, used to build the background gene
#'            list. Must have 'Protein' or 'Protein_description'.
#' @param target Character. The target namespace for gconvert (e.g., "ENSG").
#' @param numeric_ns Character. Namespace for numeric IDs.
#' @param mthreshold Numeric. g:GOSt multi-query significance threshold.
#' @param filter_na Logical. Whether gconvert should filter out NA results.
#' @param organismo Character. The organism name for gprofiler2 (e.g., "hsapiens").
#' @param custombg Logical. Whether to use a custom background (from `raw`).
#' @param platform Numeric. An indicator for the data platform (1, 2, 3, or 4).
#' @param ... Additional arguments passed to `gprofiler2::gost`.
#'
#' @return A list named `go_structures` containing:
#'         1. A `compareClusterResult` object (for plotting).
#'         2. An `enrichResult` object (for plotting).
#'         3. The raw `gost` result object (for `gostplot`).
#'
#' @export
#' @importFrom dplyr %>% bind_rows filter arrange case_when
#' @importFrom gprofiler2 gconvert gost
#' @importFrom stats p.adjust
#' @importFrom methods new
#'
#' @examples
#' # Create synthetic data
#' limma_df <- data.frame(
#'   Protein = c("P04637", "P10415", "P00533"),
#'   Protein_description = c("TP53", "CALCA", "EGFR"),
#'   expression = c("Up-regulated", "Down-regulated", "Up-regulated"),
#'   p.value = c(0.01, 0.02, 0.03)
#' )
#' raw_df <- data.frame(
#'   Protein = c("P04637", "P10415", "P00533", "P62258"),
#'   Protein_description = c("TP53", "CALCA", "EGFR", "ACTB")
#' )
#'
#' # This example requires an internet connection and will be skipped on tests
#' \donttest{
#' if (requireNamespace("gprofiler2")) {
#'   go_res <- Goterms_finder(
#'     df = limma_df,
#'     raw = raw_df,
#'     target = "ENSG",
#'     numeric_ns = "ENTREZGENE_ACC",
#'     mthreshold = 0.05,
#'     filter_na = TRUE,
#'     organismo = "hsapiens",
#'     custombg = TRUE,
#'     platform = 1
#'   )
#'   # Check the first result object
#'   if (!is.null(go_res) && inherits(go_res[[1]], "compareClusterResult")) {
#'     print(head(go_res[[1]]@compareClusterResult))
#'   }
#' }
#' }
#'

Goterms_finder <- function(df, raw,  target, numeric_ns, mthreshold, filter_na, organismo, custombg, platform, ...){
  #@ df es el dataframe que es la salida de limma.
  up <- dplyr::bind_rows(
    df %>%
      dplyr::filter(df$expression == 'Up-regulated') %>%
      dplyr::arrange(p.value)
  )
  #up <- subset(up, up$p.value <= 0.05)

  down <- dplyr::bind_rows(
    df %>%
      dplyr::filter(df$expression == 'Down-regulated') %>%
      dplyr::arrange(p.value)
  )
  #down <- subset(down, down$p.value <= 0.05)

  #Afianzamos la obtenci?n de los identificadores de nuestras prote?nas a ensemblegenome
  if (platform == 1 | platform == 2 | platform == 4){
    up_names <- gprofiler2::gconvert(up$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
    down_names <- gprofiler2::gconvert(down$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
    background <- gprofiler2::gconvert(raw$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)

  } else if (platform == 3){
    up_names <- gprofiler2::gconvert(up$Protein_description, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
    down_names <- gprofiler2::gconvert(down$Protein_description, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
    background <- gprofiler2::gconvert(raw$Protein_description, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
  }
  #Multi enrichment analysis

  if (custombg == TRUE){
    multi_gp <- gprofiler2::gost(list("up-regulated" = up_names$name, "down-regulated" = down_names$name), organism = organismo, custom_bg = background$name, ...)

  } else if (custombg == FALSE){
    multi_gp <- gprofiler2::gost(list("up-regulated" = up_names$name, "down-regulated" = down_names$name), organism = organismo, ...)
  }

  multi_gp[[1]]$adj.P.Val <- p.adjust(multi_gp[[1]]$p_value, method = "fdr") # p.adjust is base R (stats)
  gp_mod <- multi_gp$result[, c("query", "source", "term_id",    #Se extraen las columnas de inter?s
                                "term_name", "p_value","adj.P.Val", "query_size",
                                "intersection_size" , "term_size",
                                "effective_domain_size", "intersection")]
  gp_mod$ProteinRatio <- paste0(gp_mod$intersection_size, "/", gp_mod$query_size) #Creamos la columna GeneRatio

  gp_mod$BgRatio <- paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)#Creamos la columna BgRatio
  names(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.value", "adj.P.Val",
                     "query_size", "Count", "term_size", "effective_domain_size",
                     "geneID", "GeneRatio", "BgRatio")

  gp_mod$geneID <- gsub(",", "/", gp_mod$geneID)
  gp_mod$Conditions <- dplyr::case_when( gp_mod$Cluster == "up-regulated" ~ "Up-regulated",
                                         gp_mod$Cluster == "down-regulated" ~ "Down-regulated")
  gp_mod2 <- gp_mod[!duplicated(gp_mod$ID), ]
  row.names(gp_mod2) <- gp_mod2$ID

  # definimos objeto compareCluster
  # 'new' is base R, but "compareClusterResult" class is from clusterProfiler
  gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod2)
  # definimos objeto enrichResult
  # "enrichResult" class is from DOSE/enrichplot
  gp_mod_enrich <- new("enrichResult", result = gp_mod2)

  go_structures <- list(gp_mod_cluster, gp_mod_enrich, multi_gp)

  return(go_structures)
}


#' @title Create Dot Plot from Enrichment Results
#' @description Generates a dot plot from a `compareClusterResult` object,
#' faceted by condition.
#'
#' @param terms A list, where the first element (`terms[[1]]`) is a
#'              `compareClusterResult` object.
#' @param ... Additional arguments passed to `enrichplot::dotplot`.
#'
#' @return A ggplot object.
#' @export
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 facet_grid
#'
#' @examples
#' # This example requires a valid compareClusterResult object.
#' # We will create a minimal synthetic one.
#' if (requireNamespace("clusterProfiler") && requireNamespace("methods")) {
#'   dummy_data <- data.frame(
#'     Cluster = "Up-regulated", Category = "BP", ID = "GO:0006915",
#'     Description = "apoptosis", p.value = 0.01, adj.P.Val = 0.01,
#'     query_size = 50, Count = 5, term_size = 100,
#'     effective_domain_size = 10000, geneID = "P04637/P00533",
#'     GeneRatio = "5/50", BgRatio = "100/10000", Conditions = "Up-regulated"
#'   )
#'   dummy_cc_result <- methods::new("compareClusterResult",
#'                                 compareClusterResult = dummy_data)
#'   terms_list <- list(dummy_cc_result)
#'   dotplot_func(terms_list)
#' }
#'
dotplot_func <- function(terms, ...){

  enrichplot::dotplot(terms[[1]],  ...) +  ggplot2::facet_grid(.~Conditions)

}


#' @title Create g:Profiler Plot
#' @description Generates a publication-ready plot from a `gost` result object.
#'
#' @param terms A list, where the third element (`terms[[3]]`) is a
#'              `gost` result object.
#' @param ... Additional arguments passed to `gprofiler2::gostplot`.
#'
#' @return A ggplot object.
#' @export
#' @importFrom gprofiler2 gostplot
#'
#' @examples
#' # This example requires a valid gost object, which needs an internet connection
#' # See the example for Goterms_finder
#' \donttest{
#' if (requireNamespace("gprofiler2")) {
#'   # Example of how it would be called, but depends on Goterms_finder
#'   # go_res <- Goterms_finder(...)
#'   # if (!is.null(go_res)) {
#'   #   gostplot_func(go_res)
#'   # }
#'   print("This plot depends on a live query from Goterms_finder.")
#' }
#' }
#'
gostplot_func <- function(terms, ...){
  gprofiler2::gostplot(terms[[3]], ...)

}

#' @title Create Bar Plot from Enrichment Results
#' @description Generates a bar plot from an `enrichResult` object,
#' faceted by cluster.
#'
#' @param terms A list, where the second element (`terms[[2]]`) is an
#'              `enrichResult` object.
#' @param number Numeric. The number of top terms *per condition* to display.
#' @param conditions Character. A title for the plot.
#' @param ... Additional arguments passed to `enrichplot::barplot`.
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 facet_grid ylab ggtitle
#' @importFrom methods new
#'
#' @examples
#' # This example requires a valid enrichResult object.
#' # We will create a minimal synthetic one.
#' if (requireNamespace("DOSE") && requireNamespace("methods")) {
#'   dummy_data <- data.frame(
#'     Cluster = c("Up-regulated", "Down-regulated"),
#'     Category = "BP", ID = c("GO:0006915", "GO:0007049"),
#'     Description = c("apoptosis", "cell cycle"),
#'     p.value = c(0.01, 0.02), adj.P.Val = c(0.01, 0.02),
#'     query_size = c(50, 40), Count = c(5, 8), term_size = c(100, 120),
#'     effective_domain_size = 10000, geneID = c("P04637", "P10415"),
#'     GeneRatio = c("5/50", "8/40"), BgRatio = c("100/10000", "120/10000"),
#'     Conditions = c("Up-regulated", "Down-regulated")
#'   )
#'   dummy_enrich_result <- methods::new("enrichResult", result = dummy_data)
#'   terms_list <- list(NULL, dummy_enrich_result)
#'   barplot_func(terms_list, number = 1, conditions = "My Plot")
#' }
#'
#'
barplot_func <- function(terms, number, conditions, ...){


  terms[[2]]@result <- terms[[2]]@result[order(terms[[2]]@result$Count, decreasing = TRUE), ]

  #terms_display <- number * 2


  up_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Up-regulated", ][1:number,]
  down_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Down-regulated", ][1:number,]
  terms[[2]]@result <- rbind(up_regulated_subset, down_regulated_subset)


  barplot(terms[[2]], showCategory = 100, ...) +
    ggplot2::facet_grid(~Cluster) + ggplot2::ylab("Number of proteins") + ggplot2::ggtitle(conditions)

}



##################################################################################
#Interaction network analysis
#' @title Get STRING Interactions for Up-regulated Proteins
#' @description Maps up-regulated proteins to the STRING database, plots the
#' network, and returns a list of mapped proteins.
#'
#' @param df A data frame from statistical analysis (e.g., limma), must have
#'           'expression' and 'Protein' columns.
#' @param taxonid Numeric. The NCBI taxonomy ID for the species (e.g., 9606
#'                for human).
#' @param score Numeric. The minimum required interaction score (0-1000).
#'
#' @return A list containing two data frames:
#'         1. `up_mapped`: Mapped up-regulated proteins.
#'         2. `down_mapped`: Mapped down-regulated proteins.
#' @export
#' @importFrom STRINGdb STRINGdb
#' @importFrom graphics par
#' @importFrom methods new
#'
#' @examples
#' \donttest{
#' if (requireNamespace("STRINGdb")) {
#'   # Create synthetic data
#'   limma_df <- data.frame(
#'     Protein = c("P04637", "P00533"), # TP53, EGFR
#'     expression = c("Up-regulated", "Up-regulated")
#'   )
#'   # Run function
#'   interactions_list <- interactions_up(limma_df, taxonid = 9606, score = 400)
#' }
#' }
#'
interactions_up <- function(df, taxonid, score){
  #Subset de up-regulated y down-regulated
  string_db <- STRINGdb::STRINGdb$new(version="12", species=taxonid, score_threshold=score, input_directory="")

  up_regulated <- subset(df, df$expression == "Up-regulated")
  down_regulated <- subset(df, df$expression == "Down-regulated")

  up_mapped <- string_db$map(up_regulated, "Protein", removeUnmappedRows = TRUE)

  down_mapped <- string_db$map(down_regulated, "Protein", removeUnmappedRows = TRUE)

  par(mfrow=c(1,1)) # par is base R
  hits_up <- up_mapped$STRING_id[1:100]
  string_db$plot_network(hits_up)

  interactions <- list(up_mapped, down_mapped)
  return(interactions)

}


#' @title Plot STRING Interactions for Down-regulated Proteins
#' @description Maps down-regulated proteins to the STRING database and
#' plots the network.
#'
#' @param df A data frame from statistical analysis (e.g., limma), must have
#'           'expression' and 'Protein' columns.
#' @param taxonid Numeric. The NCBI taxonomy ID for the species (e.g., 9606
#'                for human).
#' @param score Numeric. The minimum required interaction score (0-1000).
#'
#' @return This function is called for its side effect of generating a plot.
#'         It does not return an object.
#' @export
#' @importFrom STRINGdb STRINGdb
#' @importFrom graphics par
#' @importFrom methods new
#'
#' @examples
#' \donttest{
#' if (requireNamespace("STRINGdb")) {
#'   # Create synthetic data
#'   limma_df <- data.frame(
#'     Protein = c("P10415", "P62258"), # CALCA, ACTB
#'     expression = c("Down-regulated", "Down-regulated")
#'   )
#'   # Run function
#'   interactions_down(limma_df, taxonid = 9606, score = 400)
#' }
#' }
#'
interactions_down <- function(df, taxonid, score){
  string_db <- STRINGdb::STRINGdb$new(version="12", species=taxonid, score_threshold=score, input_directory="")

  down_regulated <- subset(df, df$expression == "Down-regulated")

  down_mapped <- string_db$map(down_regulated, "Protein", removeUnmappedRows = TRUE)

  par(mfrow=c(1,1)) # par is base R
  hits_down <- down_mapped$STRING_id[1:100]
  string_db$plot_network(hits_down)


}

##################################################################################
#igraph analysis

#' @title Perform igraph Analysis on STRING Subnetworks
#' @description Calculates graph theory measures for up- and down-regulated
#' protein networks.
#'
#' @param interactions A list of data frames (like that returned by
#'                     `interactions_up`) containing mapped STRING IDs.
#'                     `interactions[[1]]` = up-regulated,
#'                     `interactions[[2]]` = down-regulated.
#' @param taxonid Numeric. The NCBI taxonomy ID for the species (e.g., 9606
#'                for human).
#' @param score Numeric. The minimum required interaction score (0-1000).
#'
#' @return A list of two data frames, `Upregulated` and `Downregulated`,
#'         containing graph measures for each network.
#' @export
#' @importFrom STRINGdb STRINGdb
#' @importFrom methods new
#'
#' @examples
#' \donttest{
#' if (requireNamespace("STRINGdb") && requireNamespace("igraph")) {
#'   # Create synthetic data
#'   up_mapped <- data.frame(
#'     Protein = c("TP53", "EGFR"),
#'     STRING_id = c("9606.ENSP00000269305", "9606.ENSP00000275493")
#'   )
#'   down_mapped <- data.frame(
#'     Protein = c("CALCA", "ACTB"),
#'     STRING_id = c("9606.ENSP00000302220", "9606.ENSP00000482455")
#'   )
#'   interactions_list <- list(up_mapped, down_mapped)
#'
#'   # Run analysis
#'   analysis_results <- igraph_analysis(interactions_list, 9606, 400)
#'   print(analysis_results$Upregulated)
#' }
#' }
#'
igraph_analysis <- function(interactions, taxonid, score) {
  string_db <- STRINGdb::STRINGdb$new(version = "12", species = taxonid, score_threshold = score, input_directory = "")

  hits_upregulated <- interactions[[1]]
  hits_downregulated <- interactions[[2]]

  subgraph_up_proteins <- string_db$get_subnetwork(hits_upregulated$STRING_id[1:100])
  subgraph_down_proteins <- string_db$get_subnetwork(hits_downregulated$STRING_id[1:100])

  graph_analysis_up <- calculate_graph_measures(subgraph_up_proteins, hits_upregulated)
  graph_analysis_down <- calculate_graph_measures(subgraph_down_proteins, hits_downregulated)

  list(Upregulated = graph_analysis_up, Downregulated = graph_analysis_down)
}
calculate_graph_measures <- function(subgraph, hits) {
  measures <- data.frame(
    order = igraph::vcount(subgraph),
    size = igraph::ecount(subgraph),
    density = igraph::edge_density(subgraph),
    components = igraph::count_components(subgraph),
    Clustering_coefficient = igraph::transitivity(subgraph)
  )

  degree_values <- igraph::degree(subgraph)
  top_deg_index <- order(degree_values, decreasing = TRUE)[1:10]
  top_deg_proteins <- hits$Protein[hits$STRING_id %in% igraph::V(subgraph)$name[top_deg_index]]

  betweenness_values <- igraph::betweenness(subgraph, directed = TRUE, weights = NA)
  top_betweenness_index <- order(betweenness_values, decreasing = TRUE)[1:10]
  top_betweenness_proteins <- hits$Protein[hits$STRING_id %in% igraph::V(subgraph)$name[top_betweenness_index]]

  eigen_values <- igraph::eigen_centrality(subgraph, directed = TRUE, weights = NA)$vector
  top_eigen_index <- order(eigen_values, decreasing = TRUE)[1:10]
  top_eigen_proteins <- hits$Protein[hits$STRING_id %in% igraph::V(subgraph)$name[top_eigen_index]]

  closeness_values <- igraph::closeness(subgraph, mode = "all")
  closeness_values[is.infinite(closeness_values)] <- 0
  top_closeness_index <- order(closeness_values, decreasing = TRUE)[1:10]
  top_closeness_proteins <- hits$Protein[hits$STRING_id %in% igraph::V(subgraph)$name[top_closeness_index]]

  graph_analysis <- data.frame(
    measures,
    top_deg_proteins = paste(top_deg_proteins, collapse = ";"),
    top_betweenness_proteins = paste(top_betweenness_proteins, collapse = ";"),
    top_eigen_proteins = paste(top_eigen_proteins, collapse = ";"),
    top_closeness_proteins = paste(top_closeness_proteins, collapse = ";")
  )

  return(graph_analysis)
}
