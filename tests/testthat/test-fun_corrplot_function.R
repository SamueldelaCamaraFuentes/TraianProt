library(testthat)
library(TraianProt)
library(ggplot2)

# ------------------------------------------------------------------------------
# TEST: Correlation Plot
# Objective: Verify that the correlation matrix is calculated and plotted.
# ------------------------------------------------------------------------------
test_that("corrplot_function generates a ggplot with the correlation matrix", {
    df_test <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        Sample_A = c(10, 12, 14, 15),
        Sample_B = c(10.5, 11.5, 14.2, 14.9),
        Sample_C = c(2, 4, 1, 3)
    )

    metadata_test <- data.frame(
        log2_col = c("Sample_A", "Sample_B", "Sample_C"),
        group = c("Group1", "Group1", "Group2")
    )

    plot_res <- corrplot_function(df_test, metadata_test)

    expect_s3_class(plot_res, "ggplot")

    expect_equal(nrow(plot_res$data), 9)

    self_cor <- plot_res$data[plot_res$data$Var1 == "Sample_A" & plot_res$data$Var2 == "Sample_A", "value"]
    expect_equal(as.numeric(self_cor), 1)
})
