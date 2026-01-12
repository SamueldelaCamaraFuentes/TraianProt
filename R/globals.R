#' @importFrom data.table fread
#' @importFrom dplyr rowwise c_across if_else
#' @importFrom ggplot2 remove_missing
#' @importFrom pheatmap pheatmap
#' @importFrom stringr str_replace
#' @importFrom wrProteo readMaxQuantFile
#' @importFrom clusterProfiler enrichGO
NULL

utils::globalVariables(c(
  # Variables anteriores
  "group", "log2_col", "Condition", "val", "Samples",
  "p.value", "Protein", "Intensity", "Sample",
  "Potential.contaminant", "Reverse", "Only.identified.by.site",
  "Q.Value", "counts", "PC_X", "PC_Y", "Dim1", "Dim2",
  "adj.P.Val", "sca.P.Value", "sca.adj.pval", "x", "y",
  "Var1", "Var2", "value", "KEEP", "Type",
  "counts_condition1", "counts_condition2",
  "intensity_sample_name", "raw_name", "sample_name",
  "unique_peptides_col", "Count", "Group", "Intensity_scaled"
))
