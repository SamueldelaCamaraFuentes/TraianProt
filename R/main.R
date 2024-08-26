#################################### Pipeline #################################### 


#Installing CRAN packages:
if(!require(BiocManager)){install.packages("BiocManager")}
library(BiocManager)
if(!require(methods)){install.packages("methods")}
library(methods)
if(!require(colorspace)){install.packages("colorspace")}
library(colorspace)
if(!require(grid)){install.packages("grid")}
library(grid)
if(!require(utils)){install.packages("utils")}
library(utils)
if(!require(dplyr)){install.packages("dplyr")}
library(dplyr)
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
if(!require(VennDiagram)){install.packages("VennDiagram")}
library(VennDiagram)
if(!require(VIM)){install.packages("VIM")}
library(VIM)
if(!require(wrProteo)){install.packages("wrProteo")}
library(wrProteo)
if(!require(wrMisc)){install.packages("wrMisc")}
library(wrMisc)
if(!require(gplots)){install.packages("gplots")}
library(gplots)
if(!require(gprofiler2)){install.packages("gprofiler2")}
library(gprofiler2) 
if(!require(writexl)){install.packages("writexl")}
library(writexl)
if(!require(igraph)){install.packages("igraph")}
library(igraph)
if(!require(plotly)){install.packages("plotly")}
library(plotly)
if(!require(matrixStats)){install.packages("matrixStats")}
library(matrixStats)
if(!require(ggExtra)){install.packages("ggExtra")}
library(ggExtra)
if(!require(matrixStats)){install.packages("matrixStats")}
library(matrixStats)
if(!require(DEqMS)){install.packages("DEqMS")}
library(DEqMS)
if(!require(dplyr)){install.packages("dplyr")}
library(dplyr)
if(!require(corrplot)){install.packages("corrplot")}
library(corrplot)
if(!require(devtools)){install.packages("devtools")}
library(devtools)
library(diann)


#Installing bioconductor packages: 
if(!require(rsconnect)){BiocManager::install("rsconnect",update=F,ask=F)}
library(rsconnect)
if(!require(BiocGenerics)){BiocManager::install("BiocGenerics",update=F,ask=F)}
library(BiocGenerics)
if(!require(Biobase)){BiocManager::install("Biobase",update=F,ask=F)}
library(Biobase)
if(!require(S4Vectors)){BiocManager::install("S4Vectors",update=F,ask=F)}
library(S4Vectors)
if(!require(vsn)){BiocManager::install("vsn",update=F,ask=F)}
library(vsn)
if(!require(IRanges)){BiocManager::install("IRanges",update=F,ask=F)}
if(!require(AnnotationDbi)){BiocManager::install("AnnotationDbi",update=F,ask=F)}
library(AnnotationDbi)
if(!require(limma)){BiocManager::install("limma",update=F,ask=F)}
library(limma)
if(!require(qvalue)){BiocManager::install("qvalue",update=F,ask=F)}
library(qvalue)
if(!require(clusterProfiler)){BiocManager::install("clusterProfiler",update=F,ask=F)}
library(clusterProfiler)
if(!require(enrichplot)){BiocManager::install("enrichplot",update=F,ask=F)}
library(enrichplot)
if(!require(DOSE)){BiocManager::install("DOSE",update=F,ask=F)}
library(DOSE)
if(!require(STRINGdb)){BiocManager::install("STRINGdb",update=F,ask=F)}
library(STRINGdb)


packagesdir <- readline("Introduzca el directorio donde se encuentran las funciones:")
#E:/SAMUEL/Mi app/Pipeline

packagesdir <- "E:/SAMUEL/Mi app/Pipeline"

setwd(packagesdir)

source("dataprocessing.R")
source("quality_metrics.R")
source("statisticalAnalysis.R")
source("overviewfigures.R")
source("enrichment.R")
source("Network.R")

#Source dir
initialDirectory <- readline("Introduzca el directorio donde estan los resultados de MaxQuant:")
initialDirectory <- "E:/SAMUEL/Mi app/Pipeline/datasets"

#F:/Master/Trabajo Fin de Master
#E:/SAMUEL/Mi app/Pipeline/datasets

setwd(initialDirectory)


newDirectory <- readline('Introduzca el directorio donde generar los resultados:') 
#E:/SAMUEL/Mi app/Pipeline/datasets

setwd(newDirectory)

dir_create <- readline("Introduzca el nombre del directorio donde guardar los plots:")
#E:/SAMUEL/Mi app/Pipeline/datasets

dir.create(dir_create) 
finalDirectory = paste0(dir_create)
setwd(finalDirectory)

##################################################################################
#Filtrado inicial
replicas_condicion1 <- as.numeric(readline("Introduzca el numero de replicas biologicas de su condicion control:"))
if (eval(is.na(replicas_condicion1))) {
  replicas_condicion1 <- as.numeric(readline("Error message, Introduzca con numero el numero de replicas biologicas de su condicion control:"))
}
replicas_condicion2 <- as.numeric(readline("Introduzca el numero de replicas biologicas de su condicion tratamiento:"))
if (eval(is.na(replicas_condicion2))) {
  replicas_condicion2 <- as.numeric(readline("Error message, Introduzca con numero el numero de replicas biologicas de su condicion tratamiento:"))
}

input <- readline("Introduzca el tipo de intensidad a ) lfq o bien b ) int:")

cat("Decida plataforma:")
cat("1) MaxQuant.")
cat("2) MSFragger.")
cat("3) Dia-NN.")
cat("4) Proteome Discoverer")
platform <- as.numeric(readLines(con = stdin(), n = 1))
if (eval(is.na(platform))) {
  platform <- as.numeric(readline("Error message, Introduzca con numero la plataforma seleccionada:"))
}

if (platform == 1 | platform == 2){
  raw <- read.delim("proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, colClasses = "character") 
  }else if (platform == 3){
    #raw <- read.delim("Chymeric pairwaise ttets_abundanceNORMALIZED+count.txt", sep = "\t", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE) 
    raw <- read.delim("report.pg_matrix.tsv", sep = "\t", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE) 
    
  } else if (platform ==4){
    raw <- as.data.frame(readxl::read_xlsx("Chymeric pairwaise ttets_abundanceNORMALIZED+count.xlsx", sheet = "Proteins"))
  }

#Tipo de marcaje 

cat("Marcaje:")
cat("1) Label free.")
cat("2) TMT.")
labeltype <- as.numeric(readLines(con = stdin(), n = 1))

cat("Organismo:")
cat("1) Candida albicans.")
cat("2) Otro.")
organism <- as.numeric(readLines(con = stdin(), n = 1))

if (eval(is.na(platform))) {
  organism <- as.numeric(readline("Error message, Introduzca con numero la opcion de organismo:"))
}


first_condition <- readline("Introduzca la expresion regular de su primera condicion:") #WT[0-9]
second_condition <- readline("Introduzca la expresion regular de su segunda condicion:") #WT_H2O2
df<- quick_filtering(raw, input, platform, organism, first_condition, second_condition)

#We set as datafrane in order to avoid problem of merging tables onwards. 
if (platform == 3){
  df <- as.data.frame(df)
}

##################################################################################
#Obtencion Log2 names

LOG2.names <- obtain_LOG.names(df)

##################################################################################

conditions1 <- c(first_condition, second_condition)
condition1_names <- grep(conditions1[1], LOG2.names, value = TRUE)
condition2_names <- grep(conditions1[2], LOG2.names, value = TRUE)
cond.names <- c(condition1_names,condition2_names)

#Exclusive proteins

if (labeltype == 2){  #We add this detail in order to get rid of the proteins that contain NAs in all the measurements with TMT label
  df <- df[rowSums(is.na(df[, cond.names])) != length(cond.names), ]
}

unique_proteins <- obtain_unique_proteins(df, conditions1, LOG2.names, replicas_condicion1, replicas_condicion2)




##################################################################################
#Venn Diagram

venn <- venn_diagram(df, unique_proteins, label1 = "Treatment", label2 = "Control", color1 = "blue", color2 = "maroon")
grid.draw(venn)

##################################################################################
#Proteins identified

identify_proteins(df, cond.names, platform, replicas_condicion1, replicas_condicion2) 
##################################################################################
#Preprocesamiento

min_number <- as.numeric(readline("Introduzca el minimo numero de replicas biologicas en los que encontrar cada proteina:"))
if (eval(is.na(min_number))) {
  min_number <- as.numeric(readline("Error message, Introduzca el minimo numero de replicas biologicas en los que encontrar cada proteina:"))
} else if (min_number > replicas_condicion1){
  min_number <- as.numeric(readline("Error message, Introduzca un numero que no sea mayor al numero de replicas:"))
} else if (min_number > replicas_condicion2){
  min_number <- as.numeric(readline("Error message, Introduzca un numero que no sea mayor al numero de replicas:"))
}
min_count <- c(min_number,min_number)

#Filtrado

df.F <- filter_valids(df, unique_proteins, conditions1, min_count, at_least_one <- TRUE, LOG2.names, labeltype)


#Normalizacion e imputacion

cat("Decida los pasos a proceder:")
cat("1) Normalizar e imputar.")
cat("2) Solo imputar.")
cat("3) Solo normalizar")
cat("4) No normalizar y no imputar")
choice <- as.numeric(readLines(con = stdin(), n = 1))

if (eval(is.na(choice))) {
  choice <- as.numeric(readline("Error message, Introduzca un numero:"))
}
if (choice == 1){
  df.F <- median_centering(df.F, LOG2.names)
  cat("Decida un metodo para imputar:")
  cat("1) Imputacion distribucion normal.")
  cat("2) Imputacion por metodo K-Nearest Neighbours.")
  imp_choice <- as.numeric(readLines(con = stdin(), n = 1))
  if (eval(is.na(imp_choice))) {
    imp_choice <- as.numeric(readline("Error message, Introduzca un numero:"))
  }
  if (imp_choice == 1){
    df.FNI <- impute_data(df.F, LOG2.names)
  } else if (imp_choice == 2){
    df.FNI <- impute_KNN_data(df.F, LOG2.names, k = 5)
  }
} else if (choice == 2){
  cat("Decida un metodo para imputar:")
  cat("1) Imputacion distribucion normal.")
  cat("2) Imputacion por metodo K-Nearest Neighbours.")
  imp_choice2 <- as.numeric(readLines(con = stdin(), n = 1))
  if (imp_choice2 == 1){
    df.FNI <- impute_data(df.F, LOG2.names)
  } else if (imp_choice2 == 2){
    df.FNI <- impute_KNN_data(df.F, LOG2.names, k = 5)
  }
  
} else if (choice ==3){
  df.FNI <- median_centering(df.F, LOG2.names)
} else if (choice ==4){
  df.FNI <- df.F
}

total <- bind_rows(df.FNI, as.data.frame(unique_proteins[1], check.names = FALSE))
total_dataset <- bind_rows(total, as.data.frame(unique_proteins[2], check.names = FALSE))

# Quality metrics  

#CV2
plotCV2(df.FNI[,LOG2.names],  trend = TRUE, main = "Dispersion check", cex = 0.2, pch = 16, xlab="Average log-intensity", ylab=expression("Relative standard deviation"))

#Boxplot
filename = "Boxplot.tiff"
tiff(filename = filename, units="in", width=12, height=10, res=400) #, width = 12, height = 10, units = "in", res = 400
boxplot <- boxplot_function(df.FNI, cond.names, cex.axis = 0.5)
dev.off()
boxplot <- boxplot_function(df.FNI, cond.names, cex.axis = 0.5)

#Check NAs
imputation_state(df.F, df.FNI, cond.names, cex = 1)

#Historgram
histogram(df.FNI, cond.names, color = "blue")

#Scatter plot
corr_plot <- scatterplot_function(df.FNI, cond.names[1], cond.names[2])

filename1 <- "Scatterplot_definitivo.tiff"
tiff(filename = filename1, units="in", width=12, height=10, res=400)
ggMarginal(corr_plot, type = "densigram")
dev.off()

#QQ plot
filename2 <- "qqplot.tiff"
tiff(filename = filename2, units="in", width=12, height=10, res=400)
qqplot_function(df.FNI, cond.names[2], cond.names[3], color = "blue")
dev.off()

#Correlation plot 
filename3 <- "correlation_plot.tiff"
tiff(filename = filename3, units="in", width=12, height=10, res=400)
corrplot_function(df.FNI[LOG2.names], "shade" )
dev.off()


##################################################################################
#Expresion diferencial

cat("Test estadistico:")
cat("1) Simple T-test.")
cat("2) Limma.")
test <- as.numeric(readLines(con = stdin(), n = 1))

cat("Control por PSMs, escriba TRUE o FALSE:")
cat("1) TRUE.")
cat("2) FALSE.")
psms <- as.logical(readLines(con = stdin(), n = 1))

cat("Pareado, escriba TRUE o FALSE:")
cat("1) TRUE.")
cat("2) FALSE.")
paired <- as.logical(readLines(con = stdin(), n = 1))

cat("Tipo:")
cat("1) P value.")
cat("2) Adj p value.")
statval <- as.numeric(readLines(con = stdin(), n = 1))

cat("Tipo:")
cat("1) Solo comunes")
cat("2) Incluir todo o nada.")
way <- as.numeric(readLines(con = stdin(), n = 1))

logfcup <- as.numeric(readline("Cota superior de logFC:"))
logfcdown <- as.numeric(readline("Cota inferior de logFC:"))
sig <- as.numeric(readline("Corte de p value:"))

limma <- statistical_analysis(df.FNI, LOG2.names, test =test, paired = paired, replicas_condicion1, replicas_condicion2, first_condition, second_condition, logfcup, logfcdown, sig, adjval = "fdr", statval =  statval, unique_proteins = unique_proteins, way = way, psms = psms, platform = platform) #LOG2.WT1, LOG2.WT_H2O2_1, LOG2.Prn1_1 LOG2.Prn1_H2O2_1

limma <- limma[-6]
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


write_xlsx(limma, path = "AD_Colon_Tratamiento_vs_Control.xlsx")
##################################################################################
#Representaciones
conditions <- readline("Introduzca el titulo de sus condiciones a enfrentar:") #"prn1\u2206 vs prn1\u2206 H2O2 10mM"

cat("Grosor de puntos:")
label <- as.numeric(readLines(con = stdin(), n = 1))

cat("Parametro:")
cat("1) P value.")
cat("2) Q value.")
param <- as.numeric(readLines(con = stdin(), n = 1))


volcano <- volcano_plot(limma, conditions, label, param, psms)
print(volcano)

my_pca <- pca(df.FNI, cond.names, replicas_condicion1, replicas_condicion2)

heatmap_plot <- my_heatmap(limma, cond.names, conditions)

differential_heatmap_plot <- my_heatmap_differential(limma, df.FNI, cond.names, conditions)

file = "Heatmap.tiff"
tiff(file, width = 12, height = 10, units = "in", res = 400)
my_heatmap_differential(limma, df.FNI, cond.names, conditions)
dev.off()

Diferential_boxplot(df.FNI, first_condition = first_condition, second_condition = second_condition, protein = "POR1", LOG2.names)


##################################################################################
#Analisis funcional

id_target <- readline("Introducir el tipo de identificador proteico al que mapear:") #ENSG
organismo <- readline("Introducir organismo:") #calbicans
threshold <- as.numeric(readline("Introduzca un valor umbral para la significancia:"))
if (eval(is.na(threshold))) {
  threshold <- as.numeric(readline("Error message, Introduzca un valor umbral numerico para la significancia:"))
}
Go_terms <- Goterms_finder(limma, target = id_target, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organism = organismo, user_threshold = threshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC", "CORUM"))

#Go_terms <- ggallus_enrichment(limma, target = id_target, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organism = organismo, user_threshold = threshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC", "CORUM"))

write_xlsx(Go_terms[[2]]@result, path = "AF_Colon_Tratamiento_vs_Control.xlsx")

font_size <- as.numeric(readline("Introduzca el tamaño fuente:"))
if (eval(is.na(font_size))) {
  font_size <- as.numeric(readline("Error message, Introduzca un tamaño fuente numerico:"))
}
terms_number <- as.numeric(readline("Introduzca el numero de terminos a mostrar:"))
if (eval(is.na(terms_number))) {
  terms_number <- as.numeric(readline("Error message, Introduzca un numero de terminos a mostrar numericos:"))
}

filename = "Dotplot.tiff"
tiff(filename = filename, units="in", width=12, height=10, res=400) #, width = 12, height = 10, units = "in", res = 400
dotplot_func(Go_terms, x = "GeneRatio", title = conditions, split = "Conditions", font.size = font_size, showCategory = terms_number, color = "adj.P.Val")
dev.off()

p <- gostplot_func(Go_terms, interactive = FALSE)

#Descarga del grafico de manhattan
publish_gostplot(p, filename = "manhattan.tiff")

#Descarga de la tabla con los terminos
publish_gosttable(Go_terms[[3]]$result, filename = "goterms.tiff")

filename = "Barplot.tiff"
tiff(filename = filename, units="in", width=12, height=10, res=400) #, width = 12, height = 10, units = "in", res = 400
barplot_func(Go_terms, number = terms_number, conditions, font.size = 14)
dev.off()

##################################################################################
#Analisis de interaccion de proteinas

taxonid <- as.numeric(readline("Introduzca el taxonid:"))
score <- as.numeric(readline("Introduzca puntuacion:"))

string_db <- STRINGdb$new(version="11.5", species=taxonid, score_threshold=score, input_directory="", protocol="http")

interactions_analysis <- interactions(limma)

graph_analysis <- igraph_analysis(interactions_analysis, taxonid, score)
