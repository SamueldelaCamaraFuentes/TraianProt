library(testthat)
library(TraianProt)

# ------------------------------------------------------------------------------
# TEST: Median Centering
# Objective: Verify that columns are correctly centered to median = 0,
# and that non-finite values (like Inf) are properly converted to NA.
# ------------------------------------------------------------------------------
test_that("median_centering correctly centers values and handles Inf/NA", {
    df_math <- data.frame(
        Protein = c("P1", "P2", "P3"),
        S1 = c(10, 20, 30),
        S2 = c(1, Inf, 3)
    )

    log_cols <- c("S1", "S2")

    res_centered <- median_centering(df = df_math, LOG2.names = log_cols)

    expect_equal(res_centered$S1, c(-10, 0, 10))

    expect_true(is.na(res_centered$S2[2]))

    expect_equal(res_centered$S2[1], -1)
    expect_equal(res_centered$S2[3], 1)

    expect_equal(median(res_centered$S1, na.rm = TRUE), 0)
    expect_equal(median(res_centered$S2, na.rm = TRUE), 0)
})
