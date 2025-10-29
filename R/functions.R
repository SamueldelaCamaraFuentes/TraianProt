
#################################### helper functions ####################################


##################################################################################
#Quick filtering

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

preimputation_state <- function(df, cond.names){

  par(mfrow=c(2,2)) # par is base R
  VIM::aggr(df[, cond.names], delimiter = "NA", labels = names(df[cond.names]), cex.axis = 0.4)
}

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

histogram <- function(df, colname, color, title){
  hist(as.numeric(df[,colname]), main = title, xlab = 'Intensity-value', col = color)


}

scatterplot_function <- function(df, colname1, colname2){
  corr_plot <- ggplot2::ggplot(df, ggplot2::aes(df[,colname1], df[,colname2])) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = colname2, x = colname1)

  ggExtra::ggMarginal(corr_plot, type = "densigram")

}

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

# Make sure you have these libraries loaded in your app
# library(reshape2) # for the melt function
# library(ggplot2)

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
#Limma function
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
pca <- function(x, metadata, selected_conditions){

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
  tdata<-t(x[group])

  pc<-prcomp(tdata)
  eigenValues <- pc$sdev ^2
  totalEigenValue <- sum(eigenValues)
  pcContribution <- 100 * eigenValues / totalEigenValue
  contribution1 <- round(100 * eigenValues[1] / totalEigenValue)
  contribution2 <- round(100 * eigenValues[2] / totalEigenValue)
  my_colors <-  rep(c("salmon", "blue"), c(rep1, rep2))

  plot(pc$x[,1], pc$x[,2], pch=16, cex = 3,  col = my_colors, main="Principal Component Analysis", xlab=paste0("PC1 [", contribution1, "%]"), ylab=paste0("PC2 [", contribution2, "%]"))
  text(pc$x[,1], pc$x[,2],labels=rownames(pc$x),pos=3,offset=0.4,cex=0.7)
}

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

my_heatmap <- function(data, cond.names, title){

  heatmap(as.matrix(data[cond.names]), labRow = data$Protein, main = title, Rowv = NULL, Colv = NA, col = gplots::greenred(75), cexCol = 0.6)



}
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

dotplot_func <- function(terms, ...){

  enrichplot::dotplot(terms[[1]],  ...) +  ggplot2::facet_grid(.~Conditions)

}

gostplot_func <- function(terms, ...){
  gprofiler2::gostplot(terms[[3]], ...)

}

barplot_func <- function(terms, number, conditions, ...){


  terms[[2]]@result <- terms[[2]]@result[order(terms[[2]]@result$Count, decreasing = TRUE), ]

  #terms_display <- number * 2


  up_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Up-regulated", ][1:number,]
  down_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Down-regulated", ][1:number,]
  terms[[2]]@result <- rbind(up_regulated_subset, down_regulated_subset)


  barplot(terms[[2]], showCategory = 100, ...) +
    ggplot2::facet_grid(~Cluster) + ggplot2::ylab("Number of proteins") + ggplot2::ggtitle(conditions)

}

barplot_func <- function(terms, number, conditions, ...){


  terms[[2]]@result <- terms[[2]]@result[order(terms[[2]]@result$Count, decreasing = TRUE), ]

  #terms_display <- number * 2


  up_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Up-regulated", ][1:number,]
  down_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Down-regulated", ][1:number,]
  terms[[2]]@result <- rbind(up_regulated_subset, down_regulated_subset)


  barplot(terms[[2]], showCategory = 100, ...) +
    ggplot2::facet_grid(~Cluster) + ggplot2::ylab("Number of proteins") +ggtitle(conditions)

}


##################################################################################
#Interaction network analysis

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
