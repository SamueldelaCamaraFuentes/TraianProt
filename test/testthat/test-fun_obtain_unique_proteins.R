library(testthat)
library(TraianProt)

# ------------------------------------------------------------------------------
# TEST 2: Unique Proteins Identification
# Objective: Verify that proteins present in one group but missing (Inf) in the other are identified as unique.
# ------------------------------------------------------------------------------
test_that("obtain_unique_proteins correctly identifies unique proteins", {
  
  # 1. Mock Data Setup
  # Prot1: Present in Group A (value 10), absent in Group B (Inf/Missing) -> Unique to A
  # Prot2: Present in both groups -> Not unique
  df_test <- data.frame(
    Protein = c("Prot1", "Prot2"),
    # Group A (2 replicates)
    LOG2.A1 = c(10, 10),    
    LOG2.A2 = c(10, 10),
    # Group B (2 replicates) - Inf means missing
    LOG2.B1 = c(Inf, 10),    
    LOG2.B2 = c(Inf, 10)
  )
  
  # Mock Metadata (Must include the new replicates)
  meta_test <- data.frame(
    log2_col = c("LOG2.A1", "LOG2.A2", "LOG2.B1", "LOG2.B2"),
    group = c("GroupA", "GroupA", "GroupB", "GroupB")
  )
  
  # 2. Run the function
  resultado <- obtain_unique_proteins(
    df = df_test, 
    metadata = meta_test, 
    selected_conditions = c("GroupB", "GroupA")
  )
  
  # 3. Assertions
  unique_A <- resultado[[1]]
  
  expect_true("Prot1" %in% unique_A$Protein) 
  expect_false("Prot2" %in% unique_A$Protein) 
})
