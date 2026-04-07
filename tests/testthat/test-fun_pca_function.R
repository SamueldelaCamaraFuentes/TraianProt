library(testthat)
library(TraianProt)
library(ggplot2)

# ------------------------------------------------------------------------------
# TEST: PCA Plot Generation
# Objective: Verify that the pca function correctly computes principal components,
# returns a valid ggplot object, and handles out-of-bounds PC requests.
# ------------------------------------------------------------------------------
test_that("pca function returns a valid ggplot and handles errors", {
    df_pca <- data.frame(
        CondA_1 = c(10, 12, 11, 15, 14),
        CondA_2 = c(10.5, 12.1, 11.2, 14.8, 14.1),
        CondB_1 = c(20, 22, 21, 25, 24),
        CondB_2 = c(20.5, 22.1, 21.2, 24.8, 24.1)
    )
    rownames(df_pca) <- paste0("Prot", 1:5)

    meta_pca <- data.frame(
        log2_col = c("CondA_1", "CondA_2", "CondB_1", "CondB_2"),
        group = c("Group_A", "Group_A", "Group_B", "Group_B")
    )

    conditions <- c("Group_B", "Group_A")

    res_plot <- pca(
        x = df_pca,
        metadata = meta_pca,
        selected_conditions = conditions,
        pc_x = 1,
        pc_y = 2
    )

    expect_s3_class(res_plot, "ggplot")
    expect_equal(nrow(res_plot$data), 4)

    expect_error(
        pca(df_pca, meta_pca, conditions, pc_x = 10, pc_y = 2),
        "Selected principal component index is out of bounds."
    )
})
