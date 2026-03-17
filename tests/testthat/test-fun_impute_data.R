library(testthat)
library(TraianProt)

# ------------------------------------------------------------------------------
# TEST: Imputation of missing values
# Objective: Verify that impute_data replaces NA/Inf, generates the tracking
# 'impute' columns, and uses the correct logical tracking.
# ------------------------------------------------------------------------------
test_that("impute_data handles NAs, creates tracking columns and imputes values", {
    df_missing <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        LOG2.S1 = c(20.5, NA, 19.8, 22.1), # P2 tiene un NA
        LOG2.S2 = c(21.0, 18.5, Inf, 22.0), # P3 tiene un Inf
        KEEP = c(TRUE, TRUE, TRUE, TRUE) # CRÍTICO: Tu función exige la columna KEEP
    )

    log_names <- c("LOG2.S1", "LOG2.S2")

    df_imputed <- impute_data(df = df_missing, LOG2.names = log_names)


    expect_true("impute.S1" %in% colnames(df_imputed))
    expect_true("impute.S2" %in% colnames(df_imputed))


    expect_true(df_imputed$impute.S1[2])

    expect_false(df_imputed$impute.S1[1])

    expect_true(all(is.finite(df_imputed$LOG2.S1)))
    expect_true(all(is.finite(df_imputed$LOG2.S2)))


    expect_equal(df_imputed$LOG2.S1[1], 20.5)
})
