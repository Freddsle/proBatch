test_that("constructor (wide) builds a valid ProBatchFeatures object", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    # basic construct
    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    expect_s4_class(pbf, "ProBatchFeatures")
    expect_true(methods::is(pbf, "QFeatures"))
    expect_true(validObject(pbf))

    # initial state
    log <- get_operation_log(pbf)
    expect_s4_class(log, "DataFrame")
    expect_identical(nrow(log), 0L)
    expect_true(is.list(log$params)) # list column
    expect_true(inherits(log$timestamp, "POSIXct"))
    expect_identical(get_chain(pbf), character())
    expect_identical(get_chain(pbf, as_string = TRUE), "")
    expect_identical(pb_pipeline_name(pbf), "raw")

    # assay accessors
    expect_identical(pb_current_assay(pbf), "feature::raw")
    m <- pb_assay_matrix(pbf) # default "intensity"
    expect_true(is.matrix(m))
    expect_equal(dim(m), dim(example_proteome_matrix))
    expect_equal(rownames(m), rownames(example_proteome_matrix))
    expect_equal(colnames(m), colnames(example_proteome_matrix))

    # colData alignment: rownames should match sample order
    se <- pbf[["feature::raw"]]
    expect_equal(rownames(SummarizedExperiment::colData(se)), colnames(example_proteome_matrix))
})

test_that("constructor errors on malformed inputs (wide)", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    # no colnames => error
    dm <- example_proteome_matrix
    colnames(dm) <- NULL
    expect_error(
        ProBatchFeatures(dm, example_sample_annotation, sample_id_col = "FullRunName"),
        "must have column names"
    )

    # sample_id_col not present => error
    expect_error(
        ProBatchFeatures(example_proteome_matrix, example_sample_annotation, sample_id_col = "NOPE"),
        "not found in sample_annotation"
    )

    # missing sample in annotation => error
    ann2 <- example_sample_annotation
    # drop one sample row
    ann2 <- ann2[ann2$FullRunName != colnames(example_proteome_matrix)[1], , drop = FALSE]
    expect_error(
        ProBatchFeatures(example_proteome_matrix, ann2, sample_id_col = "FullRunName"),
        "Sample annotation missing for"
    )

    # NA/empty IDs => error
    ann3 <- example_sample_annotation
    ann3$FullRunName[1] <- NA
    expect_error(
        ProBatchFeatures(example_proteome_matrix, ann3, sample_id_col = "FullRunName"),
        "NA/empty"
    )
})

test_that("constructor (long) delegates to long_to_matrix", {
    data(example_proteome, package = "proBatch") # LONG
    data(example_sample_annotation, package = "proBatch") # SA

    # Build from LONG
    pbf_long <- ProBatchFeatures_from_long(
        df_long = example_proteome,
        sample_annotation = example_sample_annotation,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity",
        name = "raw"
    )

    expect_s4_class(pbf_long, "ProBatchFeatures")
    expect_true(validObject(pbf_long))
    expect_identical(pb_current_assay(pbf_long), "feature::raw")
})

test_that("pb_as_long reuses matrix_to_long and round-trips vs direct call", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    # via helper
    long_a <- pb_as_long(
        pbf,
        feature_id_col = "peptide_group_label",
        sample_id_col  = "FullRunName",
        measure_col    = "Intensity"
    )

    # direct call (reference)
    se <- pbf[["feature::raw"]]
    long_b <- proBatch::matrix_to_long(
        data_matrix       = pb_as_wide(pbf),
        sample_annotation = as.data.frame(SummarizedExperiment::colData(se)),
        feature_id_col    = "peptide_group_label",
        measure_col       = "Intensity",
        sample_id_col     = "FullRunName"
    )

    # Compare after sorting; we only check essential columns
    keep <- c("peptide_group_label", "FullRunName", "Intensity")
    long_a2 <- long_a[, keep]
    long_b2 <- long_b[, keep]

    ord_a <- do.call(order, long_a2[, c("peptide_group_label", "FullRunName")])
    ord_b <- do.call(order, long_b2[, c("peptide_group_label", "FullRunName")])

    expect_equal(long_a2[ord_a, , drop = FALSE], long_b2[ord_b, , drop = FALSE])
})

test_that("internal logging helper updates oplog and chain", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    pbf <- ProBatchFeatures(example_proteome_matrix, example_sample_annotation, "FullRunName", "raw")

    # use internal helper via triple colon
    pbf <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "log",
        fun = "log_transform_dm",
        from = "feature::raw",
        to = "feature::raw_log",
        params = list(base = 2),
        pkg = "proBatch"
    )

    expect_true(validObject(pbf))
    log <- get_operation_log(pbf)
    expect_identical(nrow(log), 1L)
    expect_identical(as.character(log$step), "log")
    expect_identical(get_chain(pbf), "log")
    expect_identical(pb_pipeline_name(pbf), "raw")

    # add another step
    pbf <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "medianNorm",
        fun = "normalize_data_dm",
        from = "raw_log",
        to = "raw_log_mednorm",
        params = list(method = "medianCentering"),
        pkg = "proBatch"
    )
    expect_identical(get_chain(pbf), c("log", "medianNorm"))
    expect_identical(pb_pipeline_name(pbf), "raw")
})

test_that("internal add-assay-with-link links one-to-one when rows match", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    pbf <- ProBatchFeatures(example_proteome_matrix, example_sample_annotation, "FullRunName", name = "raw")

    # create a new assay with identical rownames -> should link 1:1
    se_new <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = pb_as_wide(pbf)),
        colData = SummarizedExperiment::colData(pbf[["feature::raw"]])
    )

    pbf2 <- proBatch:::.pb_add_assay_with_link(pbf, se = se_new, to = "feature::raw_copy", from = "feature::raw")
    expect_true(validObject(pbf2))
    expect_true(isTRUE(S4Vectors::metadata(pbf2)$linked_last))

    # Now modify rownames to break 1:1 -> should not link
    se_mod <- se_new
    rn <- rownames(se_mod)
    rn[1] <- paste0(rn[1], "_X")
    rownames(se_mod) <- rn

    pbf3 <- proBatch:::.pb_add_assay_with_link(pbf2, se = se_mod, to = "feature::raw_changed", from = "feature::raw_copy")
    expect_false(isTRUE(S4Vectors::metadata(pbf3)$linked_last))
})

test_that("show() prints chain and step count", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    pbf <- ProBatchFeatures(
        example_proteome_matrix,
        example_sample_annotation,
        "FullRunName",
        name = "raw"
    )

    out0 <- paste(capture.output(show(pbf)), collapse = "\n")
    expect_match(out0, "Processing chain: unprocessed data")

    # add a dummy log entry so show() reflects it
    pbf <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "log",
        fun = "log_transform_dm",
        from = "raw",
        to = "raw_log",
        params = list(base = 2),
        pkg = "proBatch"
    )

    out1 <- paste(capture.output(show(pbf)), collapse = "\n")
    expect_match(out1, "Processing chain: log")
    expect_match(out1, "Steps logged: 1")
})

test_that("validity checks fail on malformed oplog", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    pbf <- ProBatchFeatures(example_proteome_matrix, example_sample_annotation, "FullRunName", "raw")

    # Remove a required column
    bad <- get_operation_log(pbf)
    bad$fun <- NULL
    pbf_bad <- pbf
    pbf_bad@oplog <- bad
    expect_error(validObject(pbf_bad), "oplog lacks required columns")

    # Wrong types for params/timestamp
    bad2 <- S4Vectors::DataFrame(
        step      = "x",
        fun       = "f",
        from      = "a",
        to        = "b",
        params    = 1L, # not list
        timestamp = "2020-01-01", # not POSIXct
        pkg       = "p"
    )
    pbf_bad2 <- pbf
    pbf_bad2@oplog <- bad2
    expect_error(validObject(pbf_bad2), "oplog\\$params|oplog\\$timestamp")
})
