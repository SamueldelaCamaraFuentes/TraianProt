library(testthat)
library(TraianProt)

# ------------------------------------------------------------------------------
# TEST: Base R Heatmap
# Objective: Verify that my_heatmap runs without errors and successfully
# returns the clustering indices.
# ------------------------------------------------------------------------------
test_that("my_heatmap returns a list with clustering indices", {
    # 1. Mock Data Setup
    df_test <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        S1 = c(10, 12, 14, 15),
        S2 = c(10.5, 11.5, 14.2, 14.9),
        S3 = c(2, 4, 1, 3)
    )

    cond_names <- c("S1", "S2", "S3")

    # 2. Run function (silently captures the list output)
    res <- my_heatmap(data = df_test, cond.names = cond_names, title = "Test")

    # 3. Assertions
    expect_type(res, "list")
    expect_true("rowInd" %in% names(res))
    expect_true("colInd" %in% names(res))
})
