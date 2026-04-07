library(testthat)
library(TraianProt)
library(ggplot2)

# ------------------------------------------------------------------------------
# TEST: Identify Proteins (Plot generation)
# Objective: Verify that the function correctly counts proteins > 0 and
# returns a valid ggplot2 object.
# ------------------------------------------------------------------------------
test_that("identify_proteins returns a valid ggplot object with correct counts", {
    metadata_test <- data.frame(
        intensity_sample_name = c("Log_C1", "Log_C2", "Log_T1", "Log_T2"),
        group = c("Control", "Control", "Treatment", "Treatment")
    )

    raw_test <- data.frame(
        Protein = c("Prot_A", "Prot_B", "Prot_C", "Prot_D"),
        Log_C1 = c(25.1, 24.2, 0, 22.0),
        Log_C2 = c(25.0, 24.5, 0, 22.2),
        Log_T1 = c(0, 26.1, 23.5, 21.9),
        Log_T2 = c(0, 26.0, 23.8, 22.1)
    )

    plot_res <- identify_proteins(
        raw = raw_test,
        metadata = metadata_test,
        platform = 1,
        selected_conditions = c("Treatment", "Control")
    )

    expect_s3_class(plot_res, "ggplot")

    expect_equal(nrow(plot_res$data), 4)

    expect_true(all(plot_res$data$Count == 3))
})
