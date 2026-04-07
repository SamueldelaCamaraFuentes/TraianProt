library(testthat)
library(TraianProt)
library(ggplot2)

# ------------------------------------------------------------------------------
# TEST: t-SNE Plot
# Objective: Verify that the t-SNE algorithm runs and generates a ggplot object
# with the correct number of samples.
# ------------------------------------------------------------------------------
test_that("tsne returns a valid ggplot object", {
    df_test <- data.frame(
        S1 = rnorm(5), S2 = rnorm(5), S3 = rnorm(5), S4 = rnorm(5),
        S5 = rnorm(5), S6 = rnorm(5), S7 = rnorm(5), S8 = rnorm(5)
    )
    rownames(df_test) <- paste0("Prot", 1:5)

    meta_test <- data.frame(
        log2_col = paste0("S", 1:8),
        group = rep(c("CondA", "CondB"), each = 4)
    )

    res_plot <- tsne(
        x = df_test,
        metadata = meta_test,
        perplexity_num = 2,
        selected_conditions = c("CondB", "CondA")
    )

    expect_s3_class(res_plot, "ggplot")

    expect_equal(nrow(res_plot$data), 8)
})
