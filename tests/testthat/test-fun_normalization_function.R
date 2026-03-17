library(testthat)
library(TraianProt)
# ------------------------------------------------------------------------------
# TEST: Normalization wrapper (wrMisc)
# Objective: Verify that normalization_func successfully applies the specified
# normalization methods without altering the dataframe structure.
# ------------------------------------------------------------------------------
test_that("normalization_func applies normalization and keeps dimensions", {
    df_raw <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        S1 = c(10, 20, 30, 40),
        S2 = c(100, 200, 300, 400)
    )

    log_cols <- c("S1", "S2")

    df_norm_median <- normalization_func(df = df_raw, LOG2.names = log_cols, method = "median")

    expect_equal(dim(df_raw), dim(df_norm_median))

    expect_false(identical(df_raw$S2, df_norm_median$S2))

    df_norm_trim <- normalization_func(df = df_raw, LOG2.names = log_cols, method = "trimMean")

    expect_equal(dim(df_raw), dim(df_norm_trim))
    expect_true(is.numeric(df_norm_trim$S1))
})
