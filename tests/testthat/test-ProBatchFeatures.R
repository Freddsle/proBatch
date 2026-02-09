test_that("constructor (wide) builds a valid ProBatchFeatures object", {
    pb_test_load_example_data()

    # basic construct
    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    expect_s4_class(pbf, "ProBatchFeatures")
    expect_true(is(pbf, "QFeatures"))
    expect_true(validObject(pbf))

    # object-level colData must preserve sample annotation (columns + values)
    expect_equal(nrow(colData(pbf)), ncol(example_proteome_matrix))
    expect_gt(ncol(colData(pbf)), 0L)
    expect_true(all(colnames(example_sample_annotation) %in% colnames(colData(pbf))))

    # Expected sample annotation in matrix-column order
    exp_sa <- example_sample_annotation[
        match(colnames(example_proteome_matrix), example_sample_annotation$FullRunName), ,
        drop = FALSE
    ]
    rownames(exp_sa) <- exp_sa$FullRunName

    # values equal (ignore attributes/factors)
    expect_equal(
        as.data.frame(colData(pbf))[, colnames(example_sample_annotation), drop = FALSE],
        exp_sa,
        ignore_attr = TRUE
    )

    # per-assay colData must equal top-level colData
    se0 <- pbf[["feature::raw"]]
    expect_equal(
        as.data.frame(colData(pbf)),
        as.data.frame(colData(se0)),
        ignore_attr = TRUE
    )

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
    expect_equal(rownames(colData(se)), colnames(example_proteome_matrix))
})

test_that("constructor errors on malformed inputs (wide)", {
    pb_test_load_example_data()

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

    # object-level colData must be non-empty and consistent with per-assay
    expect_gt(ncol(colData(pbf_long)), 0L)
    se_long <- pbf_long[["feature::raw"]]
    expect_equal(
        as.data.frame(colData(pbf_long)),
        as.data.frame(colData(se_long)),
        ignore_attr = TRUE
    )
})

test_that("as_ProBatchFeatures wraps single-assay QFeatures with renaming", {
    skip_if_not_installed("QFeatures")

    pb_test_load_example_data()

    sample_order <- match(
        colnames(example_proteome_matrix),
        example_sample_annotation$FullRunName
    )
    expect_false(anyNA(sample_order))
    cd <- S4Vectors::DataFrame(
        example_sample_annotation[sample_order, , drop = FALSE]
    )
    rownames(cd) <- example_sample_annotation$FullRunName[sample_order]

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = example_proteome_matrix),
        colData = cd
    )

    qf <- QFeatures::QFeatures(
        experiments = list(peptideRaw = se),
        colData = cd
    )

    pbf <- as_ProBatchFeatures(qf, level = "peptide")

    expect_s4_class(pbf, "ProBatchFeatures")
    expect_true(validObject(pbf))
    expect_identical(names(pbf), "peptide::raw")
    expect_identical(pb_current_assay(pbf), "peptide::raw")
    expect_identical(get_chain(pbf), character())
    expect_identical(nrow(get_operation_log(pbf)), 0L)
    expect_equal(
        pb_assay_matrix(pbf),
        as.matrix(SummarizedExperiment::assay(qf[[1]])),
        ignore_attr = TRUE
    )
    cd_pbf <- as.data.frame(colData(pbf))
    cd_qf <- as.data.frame(colData(qf))
    expect_equal(cd_pbf[names(cd_qf)], cd_qf, ignore_attr = TRUE)
    expect_true("sample_id" %in% names(cd_pbf))
    expect_identical(cd_pbf$sample_id, rownames(SummarizedExperiment::colData(qf)))

    qf2 <- QFeatures::QFeatures(
        experiments = list(peptideRaw = se),
        colData = cd
    )
    pbf2 <- as_ProBatchFeatures(qf2, level = "peptide", pipeline = "custom")
    expect_identical(names(pbf2), "peptide::custom")
})

test_that("as_ProBatchFeatures warns on multi-assay naming mismatches", {
    skip_if_not_installed("QFeatures")

    pb_test_load_example_data()

    sample_order <- match(
        colnames(example_proteome_matrix),
        example_sample_annotation$FullRunName
    )
    expect_false(anyNA(sample_order))
    cd <- S4Vectors::DataFrame(
        example_sample_annotation[sample_order, , drop = FALSE]
    )
    rownames(cd) <- example_sample_annotation$FullRunName[sample_order]

    se1 <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = example_proteome_matrix),
        colData = cd
    )
    se2 <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = example_proteome_matrix),
        colData = cd
    )

    qf_multi <- QFeatures::QFeatures(
        experiments = list(peptideRaw = se1, proteinNorm = se2),
        colData = cd
    )

    expect_warning(
        as_ProBatchFeatures(qf_multi),
        "Some assay names do not follow"
    )
})

test_that("pb_as_long reuses matrix_to_long and round-trips vs direct call", {
    pb_test_load_example_data()

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
        sample_annotation = as.data.frame(colData(se)),
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
    pb_test_load_example_data()

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

test_that("addAssay colData conflict on DateTime is avoided by using object-level colData", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    # synthetic data with POSIXct DateTime in colData
    m <- matrix(1:6, nrow = 3, dimnames = list(paste0("f", 1:3), c("S1", "S2")))
    dt <- as.POSIXct(c("2024-01-01 00:00:00", "2024-01-02 00:00:00"), tz = "UTC")
    sa <- data.frame(
        FullRunName = c("S1", "S2"), DateTime = dt, Lab = "A",
        Condition = c("X", "Y"), row.names = c("S1", "S2"), stringsAsFactors = FALSE
    )

    pbf <- ProBatchFeatures(m, sa, "FullRunName", level = "feature") # feature::raw

    # First child assay with identical colData: ok
    se1 <- SummarizedExperiment(
        assays  = list(intensity = m),
        colData = colData(pbf) # object-level colData (good)
    )
    pbf <- proBatch:::.pb_add_assay_with_link(pbf, se = se1, to = "feature::copy", from = "feature::raw")
    expect_true(isTRUE(S4Vectors::metadata(pbf)$linked_last))

    # Now create a NEW assay but (intentionally) give it assay-level colData copy and mutate DateTime
    se2 <- SummarizedExperiment(
        assays  = list(intensity = m),
        colData = colData(pbf[["feature::copy"]])
    )
    colData(se2)$DateTime[1] <- se2$DateTime[1] + 3600 # conflict

    colData(se2) <- colData(pbf)
    pbf2 <- proBatch:::.pb_add_assay_with_link(pbf, se = se2, to = "feature::copy2", from = "feature::copy")
    expect_true(validObject(pbf2))
})

test_that(".pb_add_assay_with_link links 1:1 even if row order differs", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    m <- matrix(runif(6), nrow = 3, dimnames = list(paste0("f", 1:3), c("S1", "S2")))
    sa <- data.frame(FullRunName = c("S1", "S2"), row.names = c("S1", "S2"))
    pbf <- ProBatchFeatures(m, sa, "FullRunName", name = "raw")

    m2 <- m[sample(rownames(m)), ] # shuffled rows
    se2 <- SummarizedExperiment(
        assays  = list(intensity = m2),
        colData = colData(pbf)
    )
    pbf2 <- proBatch:::.pb_add_assay_with_link(pbf, se = se2, to = "feature::copy", from = "feature::raw")
    expect_true(isTRUE(metadata(pbf2)$linked_last))
})

test_that("pb_add_level preserves object-level colData and avoids conflicts", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    m_pep <- matrix(runif(6), nrow = 3, dimnames = list(paste0("PEP", 1:3), c("S1", "S2")))
    dt <- as.POSIXct(c("2024-06-01 08:00:00", "2024-06-01 08:10:00"), tz = "UTC")
    sa <- data.frame(FullRunName = c("S1", "S2"), row.names = c("S1", "S2"))
    pbf <- ProBatchFeatures(m_pep, sa, "FullRunName", level = "peptide") # peptide::raw

    # protein matrix and mapping (many-to-one)
    m_prot <- matrix(runif(6), nrow = 3, dimnames = list(c("A", "B", "C"), c("S1", "S2")))
    map <- data.frame(Precursor.Id = rownames(m_pep), Protein.Ids = c("A", "A", "B"))
    pbf2 <- pb_add_level(
        object = pbf,
        from = "peptide::raw",
        new_matrix = m_prot,
        to_level = "protein",
        mapping_df = map,
        from_id = "Precursor.Id",
        to_id = "Protein.Ids",
        map_strategy = "first"
    )
    expect_true(validObject(pbf2))
    expect_true("protein::raw" %in% names(pbf2))
    # identical colData across assays
    expect_identical(rownames(colData(pbf2)), rownames(colData(pbf2[["protein::raw"]])))
})

test_that("pb_add_level errors if any parent (from) is missing in mapping_df", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    m_pep <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("PEP1", "PEP2"), c("S1", "S2")))
    sa <- data.frame(FullRunName = c("S1", "S2"), row.names = c("S1", "S2"))
    pbf <- ProBatchFeatures(m_pep, sa, "FullRunName", level = "peptide")

    m_prot <- matrix(1, nrow = 1, ncol = 2, dimnames = list("A", c("S1", "S2")))
    # Only PEP1 mapped; PEP2 missing
    map <- data.frame(Precursor.Id = "PEP1", Protein.Ids = "A", stringsAsFactors = FALSE)

    expect_error(
        pb_add_level(
            object = pbf,
            from = "peptide::raw",
            new_matrix = m_prot,
            to_level = "protein",
            mapping_df = map,
            from_id = "Precursor.Id",
            to_id = "Protein.Ids",
            map_strategy = "as_is"
        ),
        "Mapping incomplete|parent IDs.*no mapping"
    )
})

test_that("pb_add_level errors when selection leaves NA parents (no exact child match)", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    m_pep <- matrix(1, nrow = 1, ncol = 2, dimnames = list("PEP1", c("S1", "S2")))
    sa <- data.frame(FullRunName = c("S1", "S2"), row.names = c("S1", "S2"))
    pbf <- ProBatchFeatures(m_pep, sa, "FullRunName", level = "peptide")

    m_prot <- matrix(1, nrow = 1, ncol = 2, dimnames = list("A", c("S1", "S2")))
    # Map to non-existing child "Z" (not in r_to), so no exact match remains
    map <- data.frame(Precursor.Id = "PEP1", Protein.Ids = "Z", stringsAsFactors = FALSE)

    expect_error(
        pb_add_level(
            object = pbf,
            from = "peptide::raw",
            new_matrix = m_prot,
            to_level = "protein",
            mapping_df = map,
            from_id = "Precursor.Id",
            to_id = "Protein.Ids",
            map_strategy = "first" # or "longest" — both should fail here
        ),
        "no exact match|have no exact match|Linking failed"
    )
})

test_that("show() lists global and level-specific pipelines", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    m_pep <- matrix(
        1:6,
        nrow = 3,
        dimnames = list(paste0("p", 1:3), c("S1", "S2"))
    )
    sa <- data.frame(
        FullRunName = c("S1", "S2"),
        row.names = c("S1", "S2"),
        stringsAsFactors = FALSE
    )
    pbf <- ProBatchFeatures(m_pep, sa, "FullRunName", level = "peptide")

    m_prot <- matrix(
        c(10, 11, 20, 21),
        nrow = 2,
        dimnames = list(paste0("prot", 1:2), c("S1", "S2"))
    )
    map <- data.frame(
        Precursor.Id = c("p1", "p2", "p3"),
        Protein.Ids = c("prot1", "prot1", "prot2"),
        stringsAsFactors = FALSE
    )
    pbf <- pb_add_level(
        object = pbf,
        from = "peptide::raw",
        new_matrix = m_prot,
        to_level = "protein",
        mapping_df = map,
        from_id = "Precursor.Id",
        to_id = "Protein.Ids",
        map_strategy = "first"
    )

    pbf <- pb_transform(pbf, from = "peptide::raw", steps = "log2", store_fast_steps = TRUE)
    pbf <- pb_transform(pbf, from = "protein::raw", steps = "log2", store_fast_steps = TRUE)

    out <- paste(capture.output(show(pbf)), collapse = "\n")
    expect_match(out, "add_level\\(protein\\)_byVar")
    expect_match(out, "peptide: log2_on_raw")
    expect_match(out, "protein: log2_on_raw")
})

test_that("pb_assay_matrix and pb_as_long compute fast logged assays on demand", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(
        c(10, 20, 30, 40),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(Sample = c("s1", "s2"), row.names = c("s1", "s2"), stringsAsFactors = FALSE)
    pbf <- ProBatchFeatures(mat, sa, sample_id_col = "Sample", level = "peptide")

    pbf <- log_transform_dm(
        pbf,
        log_base = 2,
        offset = 1,
        pbf_name = "peptide::raw",
        store_fast_steps = FALSE
    )

    expect_false("peptide::log2_on_raw" %in% names(pbf))

    base_mat <- assay(pbf[["peptide::raw"]], "intensity")
    expected <- log_transform_dm.default(base_mat, log_base = 2, offset = 1)

    resolved <- pb_assay_matrix(pbf, "peptide::log2_on_raw")
    expect_equal(resolved, expected)

    long_fast <- pb_as_long(
        pbf,
        feature_id_col = "Feature",
        sample_id_col = "Sample",
        measure_col = "Intensity",
        pbf_name = "peptide::log2_on_raw"
    )
    manual_long <- matrix_to_long(
        data_matrix = expected,
        sample_annotation = as.data.frame(colData(pbf[["peptide::raw"]])),
        feature_id_col = "Feature",
        measure_col = "Intensity",
        sample_id_col = "Sample"
    )

    ord <- order(long_fast$Feature, long_fast$Sample)
    expect_equal(long_fast[ord, ], manual_long[ord, ], ignore_attr = TRUE)
})

test_that("pb_transform honors final_name without colliding with pipeline-derived assay", {
    pb_test_load_example_data()

    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    pbf <- pb_transform(
        pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    )

    pbf <- pb_transform(
        pbf,
        from = "feature::log2_on_raw",
        steps = "medianNorm",
        params_list = list(list(
            sample_annotation = example_sample_annotation,
            sample_id_col = "FullRunName"
        ))
    )

    expect_true("feature::medianNorm_on_log2_on_raw" %in% names(pbf))


    suppressWarnings(
        pbf <- pb_transform(
            pbf,
            from = "feature::log2_on_raw",
            steps = "medianNorm",
            final_name = "feature::medianNorm_zeroNA",
            params_list = list(list(
                sample_annotation = example_sample_annotation,
                sample_id_col = "FullRunName",
                fill_the_missing = 0
            ))
        )
    )

    expect_true("feature::medianNorm_on_log2_on_raw" %in% names(pbf))
    expect_true("feature::medianNorm_zeroNA" %in% names(pbf))
})
