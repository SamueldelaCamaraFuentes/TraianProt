library(testthat)
library(TraianProt)

# ------------------------------------------------------------------------------
# TEST: Error handling in power curve
# Objective: Verify that the function throws an expected error when fed with
# non-numeric columns.
# ------------------------------------------------------------------------------
test_that("traianprot_power_curve catches non-numeric input errors", {
    df_error <- data.frame(
        Protein = c("P1", "P2", "P3"),
        R1 = c("A", "B", "C"),
        R2 = c(10, 12, 11),
        R3 = c(9, 11, 10)
    )
    expect_error(
        traianprot_power_curve(
            df = df_error,
            log2_cols = c("R1", "R2", "R3"),
            foldchange = 2,
            replicatespower = 3,
            alpha_level_choice = 1,
            alpha_level = 0.05
        ),
        "Las columnas seleccionadas en 'log2_cols' deben ser numéricas."
    )
})
