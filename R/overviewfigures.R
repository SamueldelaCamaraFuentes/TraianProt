volcano_plot <- function(limma, title, label, statval, psms){
  
  #top <- readline("Introduzca cuantas proteinas quiere etiquetar:")
  if (statval == 1){
    
    logFC <- limma$logFC
    protein_ids <- limma$Protein
    
    if (psms == TRUE){
      plot <- plot_ly(data = limma, x = ~logFC, y = ~-log10(sca.P.Value), text = protein_ids,
                      type = "scatter", mode = "markers",
                      color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                      size = I(label))
    } else if( psms == FALSE){
      plot <- plot_ly(data = limma, x = ~logFC, y = ~-log10(p.value), text = protein_ids,
                      type = "scatter", mode = "markers",
                      color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                      size = I(label))
    } 
    plot <- plot %>% layout(title = title,
                            xaxis = list(title = list(text ='Log2 Fold Change')),
                            yaxis = list(title = list(text = 'Log10 p-value')))
    
    return(plot)
  } else if (statval ==2){
    
    logFC <- limma$logFC
    protein_ids <- limma$Protein
    
    if (psms == TRUE){
      plot <- plot_ly(data = limma, x = ~logFC, y = ~-log10(sca.adj.pval), text = protein_ids,
                      type = "scatter", mode = "markers",
                      color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                      size = I(label))
    } else if( psms == FALSE){
      plot <- plot_ly(data = limma, x = ~logFC, y = ~-log10(adj.P.Val), text = protein_ids,
                      type = "scatter", mode = "markers",
                      color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                      size = I(label))
    } 
    plot <- plot %>% layout(title = title,
                            xaxis = list(title = list(text ='Log2 Fold Change')),
                            yaxis = list(title = list(text = 'Log10 q-value')))
    
    return(plot)
  }
}

pca <- function(x, group, rep1, rep2){
  x <- na.omit(x[group])
  tdata<-t(x[group])
  
  pc<-prcomp(tdata)
  eigenValues <- pc$sdev ^2
  totalEigenValue <- sum(eigenValues)
  pcContribution <- 100 * eigenValues / totalEigenValue
  contribution1 <- round(100 * eigenValues[1] / totalEigenValue)
  contribution2 <- round(100 * eigenValues[2] / totalEigenValue)
  my_colors <-  rep(c("salmon", "blue"), c(rep1, rep2))

  filename <- "PCA_plot.tiff"
  tiff(filename = filename, units="in", width=12, height=10, res=400)
  plot(pc$x[,1], pc$x[,2], pch=16, cex = 5,  col = my_colors, main="Principal Component Analysis", xlab=paste0("PC1 [", contribution1, "%]"), ylab=paste0("PC2 [", contribution2, "%]"))
  text(pc$x[,1], pc$x[,2],labels=rownames(pc$x),pos=3,offset=0.4,cex=0.9)
  dev.off()
}


my_heatmap <- function(data, cond.names, title){
  
  heatmap <- heatmap(as.matrix(data[cond.names]), labRow = data$Protein, main = title, Rowv = NULL, Colv = NA, col =greenred(75), cexCol = 0.6)
  legend(x="topleft", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")

  
}

my_heatmap_differential <- function(limma, data, cond.names, title){
  top_proteins <- bind_rows(   #Nos quedamos con las prote?nas sobre e infra expresadas en las comparativas WT
    limma %>% 
      filter(limma$expression == 'Up-regulated'),  
    limma %>% 
      filter(limma$expression == 'Down-regulated') 
  )
  
  diferential <- data.frame()
  for (i in top_proteins$Protein){   #Nos quedamos con la entrada del dataframe del dataset original de cada proteina de inter?s
    new <- as.data.frame(data[data$Protein==i,])
    diferential <- rbind(diferential, new)
  }
  diferential[cond.names] <- sapply(diferential[cond.names], as.numeric)
  heatmap <- heatmap(as.matrix(diferential[cond.names]),  labRow = diferential$Protein, main =title, Rowv = NULL, Colv = NA, col =greenred(75), cexCol = 0.6)
  legend(x="left", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")
}


Diferential_boxplot <- function(df, first_condition, second_condition, protein, LOG2.names){
   
   row.names(df) <- df$Protein
   
   subset_data <- as.matrix(df[df$Protein == protein, c(LOG2.names)])
   
   cond.colums <- grep(first_condition, colnames(df[df$Protein == protein, c(LOG2.names)]))
   trat.colums <- grep(second_condition, colnames(df[df$Protein == protein, c(LOG2.names)]))
   
   
   control <- rep(first_condition, length(cond.colums))
   treatment <- rep(second_condition, length(trat.colums))
   
   
   control_samples <- colnames(subset_data)[cond.colums]
   treatment_samples <- colnames(subset_data)[trat.colums]
   
   # Combine control and treatment samples into a single vector
   all_samples <- c(control_samples, treatment_samples)
   
   
   
   # Create a dataframe containing the values for Control and Treatment
   dt <- data.frame(
      Condition = c(control, treatment),
      val = c(as.vector(subset_data[cond.colums]), as.vector(subset_data[trat.colums])),
      Samples = all_samples
   )
   
   
   ggplot(dt, aes(Condition, val))+
      theme_dose()+
      geom_jitter(aes(color = factor(Samples)),
                  size = 5, position = position_dodge(width=0.4)) +
      geom_boxplot()+
      labs(
         y = expression(log[2]~"Intensity"),
         col = "Samples") +
      scale_color_brewer(palette = "Dark2")+
      facet_wrap(~{protein}) +
   theme(axis.title.x = element_blank())
   
}
