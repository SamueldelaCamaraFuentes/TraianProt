library(testthat)
library(TraianProt)
# ------------------------------------------------------------------------------
# TEST 3: Peptide Threshold Filtering
# Objective: Verify that proteins with fewer peptides than the threshold are removed.
# ------------------------------------------------------------------------------
test_that("unique_peptides_filter respects the minimum peptide threshold", {
  
  # 1. Mock Data: 
  df_pep <- data.frame(
    Pep.S1 = c(1, 5), # Group 1 sample
    Pep.S2 = c(1, 5)  # Group 2 sample
  )
  rownames(df_pep) <- c("ProtA", "ProtB")
  
  # 2. Metadata: 
  meta_pep <- data.frame(
    unique_peptides_col = c("Pep.S1", "Pep.S2"),
    group = c("Cond1", "Cond2")
  )
  
  # 3. Unique peptides filtering
  res <- unique_peptides_filter(
    df = df_pep, 
    metadata = meta_pep, 
    number = 2,       
    min_fraction = 1  
  )
  
  # 4. Assertions
  # ProtA tiene 1 péptido -> Debe ser eliminada
  expect_false("ProtA" %in% rownames(res)) 
  
  # ProtB tiene 5 péptidos -> Debe mantenerse
  expect_true("ProtB" %in% rownames(res))  
})