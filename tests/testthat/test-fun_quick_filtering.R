
library(testthat)
library(TraianProt)
# ------------------------------------------------------------------------------
# TEST 1: Quick Data Filtering
# Objective: Verify that contaminants, reverse hits, and site-only hits are removed in MaxQuant data.
# ------------------------------------------------------------------------------
test_that("Quick_Data_Filtering correctly removes contaminants and reverse hits", {
  
  # 1. Create Mock MaxQuant Data (Synthetic Data)
  
  raw_test <- data.frame(
    Potential.contaminant = c("+", "", "", ""),   # 1st row is a contaminant
    Reverse = c("", "+", "", ""),                 # 2nd row is a reverse hit
    Only.identified.by.site = c("", "", "+", ""), # 3rd row is identified only by site
    Fasta.headers = c("P1", "P2", "P3", "orf19.4825 orf19.4825 CGDID:CAL0005336 COORDS:Ca21chr1_C_albicans_SC5314:2119345-2118881C, translated using codon table 12 (154 amino acids) Uncharacterized ORF; Mitochondrial matrix protein; required for assembly/stability of the F1 sector of mitochondria"), # 4th row is valid
    Intensity.1 = c(100, 100, 100, 1000),         # Simulated intensity column
    stringsAsFactors = FALSE
  )
  
  # Metadata required so the function can locate the intensity columns
  meta_test <- data.frame(intensity_sample_name = "Intensity.1")
  
  # 2. Run the function
  # We use 'suppressWarnings' to keep the test output clean if non-critical warnings occur
  res <- quick_filtering(
    raw = raw_test, 
    platform = 1,      # MaxQuant
    organism = 1,      # Human
    metadata = meta_test, 
    selected_conditions = NULL
  )
  
  # 3. Assertions (Verifying expectations)
  
  # Expectation: Out of 4 rows, only 1 should remain (the 4th one)
  expect_equal(nrow(res), 1) 
  
  # Verify that the "Protein" column was created and extracted correctly from the header
  expect_equal(res$Protein, "orf19.4825")
  
  # Verify that the LOG2 transformation was performed (log2(1000) is approx 9.96)
  expect_true("LOG2.1" %in% colnames(res))
  expect_equal(res$LOG2.1[1], log2(1000))
})
