##################################################################################
#Quick filtering

quick_filtering <- function(raw, input, platform, organism, condition1, condition2){
  if (platform == 1){
    
    df <- raw %>%
      filter(Potential.contaminant != "+") %>% #Nos quedamos con todos aquellos que no tengan un +
      filter(Reverse != "+") %>% #Nos quedamos con todos aquellos que no tengan un +
      filter(Only.identified.by.site != "+") #Nos quedamos con todos aquellos que no tengan un +
    
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
    if (input == "lfq"){
      intensity_names <- grep("LFQ.intensity", colnames(df), value = TRUE)
      df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en numericos
      LOG2.names <- sub("^LFQ.intensity", "LOG2", intensity_names)
      df[LOG2.names] <- log2(df[intensity_names])
    } else if (input == "int"){
      intensity_names <- grep("^Intensity", colnames(df), value = TRUE)
      df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en integers
      LOG_names <- sub("^Intensity", "LOG2", intensity_names)
      df[LOG_names] <- log2(df[intensity_names])
    }
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
  
      if (input == "lfq"){
        intensity_names <- grep("MaxLFQ.Intensity", colnames(df), value = TRUE)
        df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en numericos
        LOG2.names <- sub("MaxLFQ.Intensity", "LOG2", intensity_names)
        #df[LOG2.names] <- log2(df[intensity_names])
        df[LOG2.names] <- lapply(df[intensity_names], log2)
      } else if (input == "int"){
        intensity_names <- grep("Intensity", colnames(df), value = TRUE)
        df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en integers
        LOG2.names <- sub("Intensity", "LOG2", intensity_names)
        df[LOG2.names] <- log2(df[intensity_names])
      }
      return(df)
    } else if (platform == 3){
      df <- raw 
      colnames(df)[1] <- "Protein"
      colnames(df)[3] <- "Protein_description"
      
      intensity_columns <- grepl("\\.d", colnames(df)) #Nos quedamos con aquellas columnas que corresponden con muestras
      intensity_columns <- colnames(df[intensity_columns])
      
      #Le damos un lavado 
      new_names <- c()
      for (i in intensity_columns){
        new <- basename(i)
        new_names <- c(new_names, new)
      }
      
      colnames(df)[length(colnames(df)) - (length(intensity_columns)-1):length(colnames(intensity_columns))] <- new_names
      intensity_names <- grep("\\.d", colnames(df), value = TRUE)
      df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en integers
      LOG2.names <- sub("\\.d", ".LOG2", intensity_names)
      df[LOG2.names] <- log2(df[intensity_names])
      
      ##############################################################################################################
      DIANN_report <- diann_load("report.tsv")
      samples <- unique(DIANN_report$Run)[order(unique(DIANN_report$Run))]
      
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
      
      # Modify the names
      
      colnames(merged_unique_peptide_counts)[2:(length(LOG2.names) +1)] <- paste("Unique peptides", colnames(merged_unique_peptide_counts)[2:(length(LOG2.names) +1)], sep=" ")
      colnames(merged_unique_peptide_counts)[1] <- "Protein"
      colnames(merged_peptide_counts)[2:(length(LOG2.names) +1)] <- paste("Peptides", colnames(merged_peptide_counts)[2:(length(LOG2.names) +1)], sep=" ")
      colnames(merged_peptide_counts)[1] <- "Protein"
      
      merged_df <- merge(df, merged_unique_peptide_counts, by = "Protein", all = FALSE)
      
      df <- merge(merged_df, merged_peptide_counts, by = "Protein", all = FALSE)
      
      
      # Sample names
      conditions1 <- c(condition1, condition2)
      condition1_names <- grep(conditions1[1], samples, value = TRUE)
      condition2_names <- grep(conditions1[2], samples, value = TRUE)
      sample_names <- c(condition1_names,condition2_names) 
      
      # Identify columns for each condition
      condition_1 <- grep(first_condition, LOG2.names, value = TRUE) 
      condition_2 <- grep(second_condition, LOG2.names, value = TRUE) 
      
      # Calculate non-NA counts for each condition separately
      df <- df %>%  #Se coge la variable del codigo
        rowwise() %>%
        mutate(
          counts_condition1 = sum(!is.na(c_across(all_of(condition_1)))),
          counts_condition2 = sum(!is.na(c_across(all_of(condition_2))))
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
    } else if (platform ==4){
      df <- raw
      abundance_names <- grep("Abundance:", colnames(df), value = TRUE)
      df[abundance_names] <- sapply(df[abundance_names], as.numeric) #convertimos los valores en numericos
      LOG2.names <- sub("^Abundance:", "LOG2", abundance_names)
      df[LOG2.names] <- log2(df[abundance_names])
      
      colnames(df)[colnames(df) == "Accession"] <- "Protein"
      
      colnames(df)[colnames(df) == "Description"] <- "Protein_description"
      if (organism == 1){
        
        regex <- regexpr(".*(?=.CGDID)", df$Protein_description, perl = TRUE)
        df$Protein <- regmatches(df$Protein_description, regex)

        #Obtenemos la descripcion proteica
        regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Protein_description, perl = TRUE)
        df$Protein_description <- regmatches(df$Protein_description, regex3)
        
      }
      return(df)
    }
  

}

##################################################################################
#Obtaining LOG2 names

obtain_LOG.names <- function(df){
  
  LOG2.names <- grep("LOG2", colnames(df))
  return(colnames(df[LOG2.names]))
}

##################################################################################
#Subseting unique proteins for each condition
obtain_unique_proteins <- function(df, conditions, LOG2.names, replicas_condicion1, replicas_condicion2) {
  cond.names <- lapply(conditions, function(x) grep(x, LOG2.names, value = TRUE, perl = TRUE))
  condi_names <- unlist(cond.names)
  
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
  sufficient_presence <- infinite_sums < ceiling(replicas_condicion2 / 2)
  condition & sufficient_presence
}))
  
  df$cond1_exclusive <- cond1_exclusive
  df$cond2_exclusive <- cond2_exclusive
  
  cond1_unicas <- filter(df, df$cond1_exclusive)
  cond2_unicas <- filter(df, df$cond2_exclusive)
  
  return(list(cond1_unicas, cond2_unicas))
}

##################################################################################
#Filtering

filter_valids <- function(df, unique_proteins, conditions1, min_count, at_least_one = FALSE, LOG2.names, labeltype) {
  
  
  cond.names <- lapply(conditions1, # Sobre la lista conditions, aplicamos una funcion quedarnos con las columnas de las condiciones con datos en log
                       function(x) grep(x, LOG2.names, value = TRUE, perl = TRUE))
  print(cond.names)
  
  total_unique_proteins <- rbind(unique_proteins[[1]], unique_proteins[[2]])
  common_df <- df%>%
    filter(df$Protein != total_unique_proteins$Protein[1])
  for (i in total_unique_proteins$Protein){
    common_df <- common_df%>%
      filter(common_df$Protein != i)
  }
  
  if (labeltype == 1){
    
  
    cond.filter <- sapply(1:length(cond.names), function(i) { #En nuestro caso se itera del 1 al 2
      df2 <- common_df[cond.names[[i]]]   # Extrae las columnas de inter?s
      df2 <- as.matrix(df2)   # Lo convierte en matriz
      sums <- rowSums(is.finite(df2)) # Cuenta el numero de valores validos para cada condicion 
      sums >= min_count[i]   # Calculates whether min_count requirement is met si el n?mero de datos es igual o superior al elemento en min_count se establece true
    })
    if (at_least_one) {
      common_df$KEEP <- apply(cond.filter, 1, any)
    } else {
      common_df$KEEP <- apply(cond.filter, 1, all)
    }
    common_df <- filter(common_df, KEEP)
    
    common_df[LOG2.names] <- lapply(LOG2.names,
                                    function(x) {
                                      temp <- common_df[[x]]
                                      temp[!is.finite(temp)] = NA
                                      return(temp)
                                      
                                    })
    
    condition1_names <- cond.names[[1]] #Manejamos grupos 
    condition2_names <- cond.names[[2]]
    
    for (i in 1:nrow(common_df)) { #Probar a poner el CV multiplicado por 100 al resultado y no al denominador
      valuesA <- as.numeric(common_df[i, condition1_names])
      valuesB <- as.numeric(common_df[i, condition2_names])
      valuesC <- c(valuesA, valuesB)
      cv1 <- (sd(valuesA, na.rm = TRUE) / mean(valuesA, na.rm = TRUE)) * 100
      cv2 <- (sd(valuesB, na.rm = TRUE) / mean(valuesB, na.rm = TRUE)) * 100
      cv3 <- (sd(valuesC, na.rm = TRUE) / mean(valuesC, na.rm = TRUE)) * 100
      common_df$CV_Control[i] <- cv1
      common_df$CV_Tratamiento[i] <- cv2
      common_df$CV_TratamientoyControl[i] <- cv3
    
      }
  
  } else if (labeltype ==2){
    common_df <- df
    common_df <- remove_missing(common_df)
    
    condition1_names <- cond.names[[1]] #Manejamos grupos 
    condition2_names <- cond.names[[2]]
    
    for (i in 1:nrow(common_df)) { #Probar a poner el CV multiplicado por 100 al resultado y no al denominador
      valuesA <- as.numeric(common_df[i, condition1_names])
      valuesB <- as.numeric(common_df[i, condition2_names])
      valuesC <- c(valuesA, valuesB)
      cv1 <- (sd(valuesA, na.rm = TRUE) / mean(valuesA, na.rm = TRUE)) * 100
      cv2 <- (sd(valuesB, na.rm = TRUE) / mean(valuesB, na.rm = TRUE)) * 100
      cv3 <- (sd(valuesC, na.rm = TRUE) / mean(valuesC, na.rm = TRUE)) * 100
      common_df$CV_Control[i] <- cv1
      common_df$CV_Tratamiento[i] <- cv2
      common_df$CV_TratamientoyControl[i] <- cv3
      
    }
  }
  
  
  return(common_df)  
}

##################################################################################
#Normalization 

median_centering <- function(df, LOG2.names) {
  
  df[, LOG2.names] <- lapply(LOG2.names, #A las columnas de inter?s, se procede a calcular el LOG2-mediana de cada 
                             function(x) {
                               LOG2 <- df[[x]] #Se asigna el valor calculado
                               LOG2[!is.finite(LOG2)] <- NA   # Excluimos los missing values del c?lculo de la mediana 
                               gMedian <- median(LOG2, na.rm = TRUE)
                               LOG2 - gMedian #Se calcula el valor.
                             }
  )
  
  return(df)
}

normalization_func <- function(df, LOG2.names, method){
  if (method == "trimMean"){
    df[, LOG2.names] <- normalizeThis(dat = df[, LOG2.names], method = method, trimFa = 0.4)
  } else if (method != "trimMean"){
    df[, LOG2.names] <- normalizeThis(dat = df[, LOG2.names], method = method)
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
  set.seed(1)
  df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- df[[x]]
                             temp[!is.finite(temp)] = NA
                             
                             temp.sd <- width * sd(temp[df$KEEP], na.rm = TRUE)   # shrink sd width
                             temp.mean <- mean(temp[df$KEEP], na.rm = TRUE) - 
                               downshift * sd(temp[df$KEEP], na.rm = TRUE)   # shift mean of imputed values
                             
                             n.missing <- sum(is.na(temp))
                             temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                             return(temp)
                           })
  #write.csv(df, paste0(finalDirectory,"/Dataset_preprocesado.xlsx"))
  
  
  return(df)
}

impute_KNN_data <- function(df, LOG2.names, ...){
  
  impute.names <- sub("LOG2", "impute", LOG2.names)
  df[impute.names] <- lapply(LOG2.names, function(x) !is.finite(df[, x])) #Creamos columnas donde se ve si imputar o no
  
  #Imputaci?n
  df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- df[[x]]
                             temp[!is.finite(temp)] = NA
                             return(temp)
                             
                           })
  
  imp_knn <- kNN(df, variable = LOG2.names, ... )
  #write.csv(imp_knn, paste0(finalDirectory,"/Dataset_preprocesado_KNN.xlsx"))
  return(imp_knn)
}