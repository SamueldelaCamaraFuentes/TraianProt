library(testthat)
library(TraianProt)

# ------------------------------------------------------------------------------
# TEST: Obtain LOG2 names
# Objective: Verify that the helper function correctly extracts ONLY the column
# names containing the string "LOG2".
# ------------------------------------------------------------------------------
test_that("obtain_LOG.names correctly identifies LOG2 columns", {
    df_test <- data.frame(
        Protein = c("ProtA", "ProtB"),
        LOG2.Control_1 = c(22.1, 23.4),
        Intensity.Control_1 = c(1000, 2000),
        LOG2.Treat_1 = c(24.1, 25.4),
        Extra_Column = c("A", "B")
    )


    res <- obtain_LOG.names(df_test)


    expect_equal(length(res), 2)
    expect_true(all(c("LOG2.Control_1", "LOG2.Treat_1") %in% res))
    expect_false("Intensity.Control_1" %in% res)
})
