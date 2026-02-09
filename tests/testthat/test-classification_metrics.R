test_that("classification metrics with provided clusters", {
    sample_ids <- paste0("S", 1:6)
    data_matrix <- matrix(
        rnorm(24),
        nrow = 4,
        dimnames = list(paste0("F", 1:4), sample_ids)
    )
    data_matrix[, 1:3] <- data_matrix[, 1:3] + 2

    sample_annotation <- data.frame(
        FullRunName = sample_ids,
        known = rep(c("A", "B"), each = 3),
        stringsAsFactors = FALSE
    )
    sample_annotation$cluster <- sample_annotation$known

    metrics <- calculate_classification_metrics(
        data_matrix,
        sample_annotation,
        known_col = "known",
        cluster_col = "cluster"
    )

    expect_s3_class(metrics, "data.frame")
    expect_equal(metrics$ARI, 1)
    expect_equal(metrics$MCC, 1)
    expect_true(is.finite(metrics$silhouette))
    expect_true(metrics$silhouette >= -1 && metrics$silhouette <= 1)

    assignments <- attr(metrics, "pb_assignments")
    expect_true(is.data.frame(assignments))
    expect_equal(assignments$known, sample_annotation$known)
    expect_equal(assignments$predicted, sample_annotation$cluster)
})

test_that("classification metrics uses kmeans fallback when no clusters supplied", {
    set.seed(1)
    sample_ids <- paste0("S", 1:6)
    data_matrix <- matrix(
        rnorm(24),
        nrow = 4,
        dimnames = list(paste0("F", 1:4), sample_ids)
    )
    data_matrix[, 1:3] <- data_matrix[, 1:3] + 1.5

    sample_annotation <- data.frame(
        FullRunName = sample_ids,
        known = rep(c("A", "B"), each = 3),
        stringsAsFactors = FALSE
    )

    metrics <- calculate_classification_metrics(
        data_matrix,
        sample_annotation,
        known_col = "known"
    )

    expect_s3_class(metrics, "data.frame")
    expect_equal(metrics$cluster_source, "kmeans")
    expect_equal(metrics$k, 2)
    expect_equal(metrics$n_clusters, 2)
    expect_true(is.finite(metrics$ARI))
    expect_true(is.finite(metrics$MCC))
    expect_true(is.finite(metrics$silhouette))
})
