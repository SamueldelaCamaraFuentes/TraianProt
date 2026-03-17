library(testthat)
library(TraianProt)
library(ggplot2)

# ------------------------------------------------------------------------------
# TEST: Differential Heatmap (ggplot)
# Objective: Verify that it plots DE proteins correctly, and handles
# gracefully the case where no DE proteins are found.
# ------------------------------------------------------------------------------
test_that("my_heatmap_differential handles DE proteins and empty results", {
    limma_de <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        expression = c("Up-regulated", "Down-regulated", "Not Significant", "Up-regulated"),
        logFC = c(2.5, -1.8, 0.5, 1.2)
    )

    df_expr <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        S1 = c(15, 5, 10, 12), S2 = c(14, 6, 11, 13),
        S3 = c(5, 15, 9, 6), S4 = c(6, 14, 10, 5)
    )

    cond_names <- c("S1", "S2", "S3", "S4")

    res_plot <- my_heatmap_differential(limma_de, df_expr, cond_names, "Test")

    expect_s3_class(res_plot, "ggplot")


    expect_equal(nrow(res_plot$data), 12)


    limma_empty <- data.frame(
        Protein = c("P1", "P2"),
        expression = c("Not Significant", "Not Significant"),
        logFC = c(0.1, -0.2)
    )

    res_empty <- my_heatmap_differential(limma_empty, df_expr, cond_names, "Empty")

    expect_s3_class(res_empty, "ggplot")

    expect_equal(nrow(res_empty$data), 0)
})
