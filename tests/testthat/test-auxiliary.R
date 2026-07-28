test_that("long_conversion_to_matrix", {
    data(example_proteome, package = "proBatch")

    rows <- c(
        which(example_proteome$peptide_group_label == "10062_NVGVSFYADKPEVTQEQK_3"),
        which(example_proteome$peptide_group_label == "10063_NVGVSFYADKPEVTQEQKK_3")
    )

    proteome <- example_proteome[rows, ]
    proteome_matrix <- long_to_matrix(proteome)

    expect_equal(
        rownames(proteome_matrix),
        c("10062_NVGVSFYADKPEVTQEQK_3", "10063_NVGVSFYADKPEVTQEQKK_3")
    )
    expect_equal(ncol(proteome_matrix), nrow(proteome) / 2)
})

test_that("long conversion rejects duplicate feature/sample keys", {
    duplicate <- data.frame(
        Feature = c("f1", "f1"),
        Sample = c("s1", "s1"),
        Intensity = c(1, 2)
    )

    expect_error(
        long_to_matrix(
            duplicate,
            feature_id_col = "Feature",
            sample_id_col = "Sample",
            measure_col = "Intensity"
        ),
        "duplicate feature/sample keys.*f1/s1"
    )
})


test_that("matrix_conversion_to_long", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    rows <- c(
        which(rownames(example_proteome_matrix) == "10062_NVGVSFYADKPEVTQEQK_3"),
        which(rownames(example_proteome_matrix) == "10063_NVGVSFYADKPEVTQEQKK_3")
    )

    peptide_matrix <- as.matrix(example_proteome_matrix[rows, ])
    proteome_long <- matrix_to_long(peptide_matrix, example_sample_annotation)

    expect_equal(
        unique(proteome_long$peptide_group_label),
        c("10062_NVGVSFYADKPEVTQEQK_3", "10063_NVGVSFYADKPEVTQEQKK_3")
    )
    expect_equal(nrow(proteome_long), ncol(peptide_matrix) * 2)

    setequal(
        colnames(proteome_long),
        c(colnames(example_sample_annotation), "peptide_group_label", "Intensity")
    )
})

test_that("matrix conversion rejects duplicate axes", {
    duplicate_features <- matrix(
        seq_len(4),
        nrow = 2,
        dimnames = list(c("f1", "f1"), c("s1", "s2"))
    )
    duplicate_samples <- matrix(
        seq_len(4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s1"))
    )

    expect_error(
        matrix_to_long(duplicate_features),
        "feature axis contains duplicate identifiers"
    )
    expect_error(
        matrix_to_long(duplicate_samples),
        "sample axis contains duplicate identifiers"
    )
})


test_that("generated_peptide_annotation", {
    data(example_proteome, package = "proBatch")

    peptide_annt <- create_peptide_annotation(example_proteome,
        feature_id_col = "peptide_group_label",
        protein_col = c("Protein")
    )

    expect_equal(colnames(peptide_annt), c("peptide_group_label", "Protein"))
})
