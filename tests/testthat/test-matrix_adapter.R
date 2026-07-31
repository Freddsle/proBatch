test_that("matrix adaptation aligns annotation and forwards method arguments", {
    input <- matrix(
        1:6,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )
    annotation <- data.frame(
        sample = c("s3", "unused", "s1", "s2"),
        shift = c(30, 99, 10, 20),
        stringsAsFactors = FALSE
    )
    observed <- new.env(parent = emptyenv())
    method <- function(data_matrix, sample_annotation, increment) {
        observed$sample_ids <- sample_annotation$sample
        observed$row_names <- rownames(sample_annotation)
        sweep(
            data_matrix,
            2,
            sample_annotation$shift + increment,
            "+"
        )
    }

    result <- pb_apply_matrix_method(
        input,
        method,
        sample_annotation = annotation,
        sample_id_col = "sample",
        increment = 1
    )

    expect_identical(observed$sample_ids, colnames(input))
    expect_identical(observed$row_names, colnames(input))
    expect_equal(result, sweep(input, 2, c(11, 21, 31), "+"))
})

test_that("matrix adaptation uses row-name fallback and registered steps", {
    input <- matrix(
        c(1, 3, 2, 6),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    annotation <- data.frame(
        batch = c("B", "A"),
        row.names = c("s2", "s1")
    )
    observed <- NULL
    identity_method <- function(data_matrix, sample_annotation) {
        observed <<- sample_annotation
        data_matrix
    }

    expect_identical(
        pb_apply_matrix_method(
            input,
            identity_method,
            sample_annotation = annotation,
            sample_id_col = "sample"
        ),
        input
    )
    expect_identical(rownames(observed), c("s1", "s2"))
    expect_identical(observed$sample, c("s1", "s2"))

    registered <- pb_apply_matrix_method(
        input,
        "medianNorm",
        sample_annotation = data.frame(
            FullRunName = colnames(input),
            row.names = colnames(input)
        )
    )
    expected <- normalize_sample_medians_dm(input)
    expect_equal(registered, expected)
    expect_equal(
        pb_apply_matrix_method(input, "log2"),
        log2(input + 1)
    )
})

test_that("matrix adaptation owns its matrix and annotation arguments", {
    input <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    method <- function(data_matrix, sample_annotation) data_matrix

    expect_error(
        pb_apply_matrix_method(
            input,
            method,
            data_matrix = input
        ),
        "adapter-owned argument.*data_matrix"
    )
    expect_error(
        pb_apply_matrix_method(
            structure(list(), class = "ProBatchFeatures"),
            method
        ),
        "must use `pb_transform\\(\\)`"
    )
})

test_that("matrix adaptation applies canonical missing-value policies", {
    input <- matrix(
        c(1, NA_real_, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    method <- function(data_matrix, sample_annotation) data_matrix

    expect_error(
        pb_apply_matrix_method(input, method),
        "Input contains missing values"
    )
    expect_error(
        pb_apply_matrix_method(
            input,
            method,
            missing = c("error", "keep", "drop_features", "fill")
        ),
        "`missing` must be one of"
    )
    for (legacy_policy in list(FALSE, 0, "remove")) {
        expect_error(
            pb_apply_matrix_method(
                input,
                method,
                missing = legacy_policy
            ),
            "`missing` must be one of"
        )
    }

    kept <- pb_apply_matrix_method(input, method, missing = "keep")
    expect_identical(kept, input)

    filled <- pb_apply_matrix_method(
        input,
        method,
        missing = "fill",
        fill_value = 0
    )
    expect_identical(unname(filled[1, 2]), 0)

    expect_warning(
        dropped <- pb_apply_matrix_method(
            input,
            method,
            missing = "drop_features"
        ),
        "removed 1 rows and 0 columns"
    )
    expect_identical(rownames(dropped), "f2")
    expect_identical(colnames(dropped), colnames(input))
})

test_that("long adaptation preserves rows and supports all retention modes", {
    input <- data.frame(
        feature = c("f2", "f1", "f2", "f1"),
        sample = c("s2", "s1", "s1", "s2"),
        intensity = c(22, 11, 21, 12),
        source = c("d", "a", "c", "b"),
        stringsAsFactors = FALSE
    )
    annotation <- data.frame(
        sample = c("s1", "s2"),
        group = c("A", "B"),
        source = c("annotation-a", "annotation-b"),
        stringsAsFactors = FALSE
    )
    subset_method <- function(data_matrix, sample_annotation) {
        data_matrix["f1", c("s2", "s1"), drop = FALSE] + 100
    }

    default <- pb_apply_matrix_method(
        input,
        subset_method,
        sample_annotation = annotation,
        feature_id_col = "feature",
        sample_id_col = "sample",
        measure_col = "intensity"
    )
    expect_identical(names(default), names(input))
    expect_identical(default$feature, c("f1", "f1"))
    expect_identical(default$sample, c("s1", "s2"))
    expect_identical(default$intensity, c(111, 112))
    expect_identical(default$source, c("a", "b"))

    all <- pb_apply_matrix_method(
        input,
        subset_method,
        sample_annotation = annotation,
        feature_id_col = "feature",
        sample_id_col = "sample",
        measure_col = "intensity",
        keep_all = "all"
    )
    expect_identical(names(all), c(names(input), "group"))
    expect_identical(all$group, c("A", "B"))
    expect_identical(all$source, c("a", "b"))

    minimal <- pb_apply_matrix_method(
        input,
        subset_method,
        sample_annotation = annotation,
        feature_id_col = "feature",
        sample_id_col = "sample",
        measure_col = "intensity",
        keep_all = "minimal"
    )
    expect_identical(
        names(minimal),
        c("feature", "sample", "intensity")
    )
    expect_identical(minimal$sample, c("s1", "s2"))
})

test_that("long adaptation rejects ambiguous keys and invalid measurements", {
    duplicate <- data.frame(
        feature = c("f1", "f1"),
        sample = c("s1", "s1"),
        intensity = c(1, 2)
    )
    character_measure <- data.frame(
        feature = "f1",
        sample = "s1",
        intensity = "one"
    )
    method <- function(data_matrix, sample_annotation) data_matrix

    expect_error(
        pb_apply_matrix_method(
            duplicate,
            method,
            feature_id_col = "feature",
            sample_id_col = "sample",
            measure_col = "intensity"
        ),
        "duplicate feature/sample keys.*f1/s1"
    )
    expect_error(
        pb_apply_matrix_method(
            character_measure,
            method,
            feature_id_col = "feature",
            sample_id_col = "sample",
            measure_col = "intensity"
        ),
        "measurement column `intensity` must be numeric",
        fixed = TRUE
    )
})

test_that("long input uses the same default missing-value policy", {
    input <- data.frame(
        feature = c("f1", "f1", "f2", "f2"),
        sample = c("s1", "s2", "s1", "s2"),
        intensity = c(1, NA_real_, 3, 4),
        token = letters[1:4],
        stringsAsFactors = FALSE
    )
    method <- function(data_matrix, sample_annotation) data_matrix

    expect_error(
        pb_apply_matrix_method(
            input,
            method,
            feature_id_col = "feature",
            sample_id_col = "sample",
            measure_col = "intensity"
        ),
        "Input contains missing values"
    )
    filled <- pb_apply_matrix_method(
        input,
        method,
        feature_id_col = "feature",
        sample_id_col = "sample",
        measure_col = "intensity",
        missing = "fill",
        fill_value = 0
    )
    expect_identical(filled$token, input$token)
    expect_equal(filled$intensity, c(1, 0, 3, 4))
})

test_that("method output must be a known ordered numeric matrix subset", {
    input <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )

    expect_error(
        pb_apply_matrix_method(
            input,
            function(data_matrix, sample_annotation) {
                rbind(data_matrix, f3 = c(5, 6))
            }
        ),
        "introduced unknown identifier"
    )
    expect_error(
        pb_apply_matrix_method(
            input,
            function(data_matrix, sample_annotation) {
                data_matrix[c("f2", "f1"), , drop = FALSE]
            }
        ),
        "preserve input feature and sample order"
    )
    expect_error(
        pb_apply_matrix_method(
            input,
            function(data_matrix, sample_annotation) {
                output <- data_matrix
                output[1, 1] <- NA_real_
                output
            }
        ),
        "Method output contains missing values"
    )
    expect_error(
        pb_apply_matrix_method(
            input,
            function(data_matrix, sample_annotation) {
                as.data.frame(data_matrix)
            }
        ),
        "Method output must be a numeric matrix",
        fixed = TRUE
    )
})

test_that("structured results preserve artifacts and adapted data", {
    matrix_input <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    matrix_method <- function(data_matrix, sample_annotation) {
        pb_step_result(
            data_matrix + 1,
            artifacts = list(model = "matrix")
        )
    }
    matrix_result <- pb_apply_matrix_method(matrix_input, matrix_method)
    expect_s3_class(matrix_result, "pb_step_result")
    expect_identical(matrix_result$data, matrix_input + 1)
    expect_identical(matrix_result$artifacts, list(model = "matrix"))

    input <- data.frame(
        feature = c("f1", "f1", "f2", "f2"),
        sample = c("s1", "s2", "s1", "s2"),
        intensity = 1:4,
        stringsAsFactors = FALSE
    )
    method <- function(data_matrix, sample_annotation) {
        pb_step_result(
            data_matrix + 10,
            artifacts = list(model = "retained")
        )
    }

    result <- pb_apply_matrix_method(
        input,
        method,
        feature_id_col = "feature",
        sample_id_col = "sample",
        measure_col = "intensity"
    )

    expect_s3_class(result, "pb_step_result")
    expect_s3_class(result$data, "data.frame")
    expect_equal(result$data$intensity, 11:14)
    expect_identical(result$artifacts, list(model = "retained"))
})
