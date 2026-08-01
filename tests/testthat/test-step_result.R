test_that("pb_step_result constructs a small validated carrier", {
    data <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    artifacts <- list(model = list(converged = TRUE), scores = c(1, 2))

    result <- pb_step_result(data, artifacts)

    expect_s3_class(result, "pb_step_result")
    expect_named(result, c("data", "artifacts"))
    expect_identical(result$data, data)
    expect_identical(result$artifacts, artifacts)
    expect_error(
        pb_step_result(data, artifacts = 1),
        "`artifacts` must be a list",
        fixed = TRUE
    )
})

test_that("step-result helpers unwrap plain and structured values", {
    data <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )

    plain <- proBatch:::.pb_step_result_parts(data)
    expect_false(plain$structured)
    expect_identical(plain$data, data)
    expect_identical(plain$artifacts, list())

    structured <- proBatch:::.pb_step_result_parts(
        pb_step_result(data, list(token = "kept"))
    )
    expect_true(structured$structured)
    expect_identical(structured$data, data)
    expect_identical(structured$artifacts, list(token = "kept"))

    malformed <- structure(
        list(data = data, artifacts = "not a list"),
        class = "pb_step_result"
    )
    expect_error(
        proBatch:::.pb_step_result_parts(malformed),
        "Invalid `pb_step_result`",
        fixed = TRUE
    )
    expect_error(
        proBatch:::.pb_step_result_parts(
            structure(list(data = data), class = "pb_step_result")
        ),
        "Invalid `pb_step_result`",
        fixed = TRUE
    )
    expect_error(
        proBatch:::.pb_step_result_parts(structure(
            list(data = data, artifacts = list(), extra = TRUE),
            class = "pb_step_result"
        )),
        "Invalid `pb_step_result`",
        fixed = TRUE
    )
})

test_that("matrix-result validation and artifact attachment are explicit", {
    data <- matrix(
        1:4,
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    parts <- proBatch:::.pb_step_result_matrix_parts(
        pb_step_result(data, list(token = TRUE)),
        "Provider result"
    )
    expect_identical(parts$data, data)

    expect_error(
        proBatch:::.pb_step_result_matrix_parts(
            as.data.frame(data),
            "Provider result"
        ),
        "Provider result must be a numeric matrix",
        fixed = TRUE
    )
    expect_error(
        proBatch:::.pb_step_result_matrix_parts(
            pb_step_result(as.data.frame(data)),
            "Provider result"
        ),
        "Provider result must contain a numeric matrix in `data`",
        fixed = TRUE
    )

    assay <- SummarizedExperiment::SummarizedExperiment(
        assays = list(data = data)
    )
    attached <- proBatch:::.pb_attach_step_artifacts(
        assay,
        list(token = "stored")
    )
    expect_identical(
        S4Vectors::metadata(attached)$pb_step_artifacts,
        list(token = "stored")
    )
    expect_error(
        proBatch:::.pb_attach_step_artifacts(assay, "invalid"),
        "`artifacts` must be a list",
        fixed = TRUE
    )
})
