library(testthat)
library(TraianProt)
library(ggplot2)

# ------------------------------------------------------------------------------
# TEST: Boxplot Generation
# Objective: Verify that boxplot_function successfully pivots the data and
# returns a valid ggplot object with boxplot geometry.
# ------------------------------------------------------------------------------
test_that("boxplot_function returns a valid ggplot boxplot", {
    df_box <- data.frame(
        Protein = c("P1", "P2", "P3", "P4"),
        CondA_1 = c(10, 11, 10.5, 12),
        CondA_2 = c(10.2, 11.1, 10.8, 11.9),
        CondB_1 = c(20, 21, 20.5, 22),
        CondB_2 = c(19.8, 21.1, 20.8, 21.9)
    )

    meta_box <- data.frame(
        log2_col = c("CondA_1", "CondA_2", "CondB_1", "CondB_2"),
        group = c("Group_A", "Group_A", "Group_B", "Group_B")
    )

    conditions <- c("Group_B", "Group_A")

    res_plot <- boxplot_function(
        df = df_box,
        metadata = meta_box,
        selected_conditions = conditions
    )

    expect_s3_class(res_plot, "ggplot")

    expect_true(inherits(res_plot$layers[[1]]$geom, "GeomBoxplot"))

    expect_equal(nrow(res_plot$data), 16)
})
