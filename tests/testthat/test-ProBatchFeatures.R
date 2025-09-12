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
        as.data.frame(SummarizedExperiment::colData(se0)),
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

    # object-level colData must be non-empty and consistent with per-assay
    expect_gt(ncol(colData(pbf_long)), 0L)
    se_long <- pbf_long[["feature::raw"]]
    expect_equal(
        as.data.frame(colData(pbf_long)),
        as.data.frame(SummarizedExperiment::colData(se_long)),
        ignore_attr = TRUE
    )
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

    # remove a DateTime column to avoid colData conflict error
    SummarizedExperiment::colData(se_mod)$DateTime <- NULL

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

test_that("pb_add_level: 1:1 linking works and uses addAssayLinkOneToOne", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    # Build a tiny peptide assay
    set.seed(1)
    m_pep <- matrix(round(rnorm(6), 3),
        nrow = 3,
        dimnames = list(paste0("pep", 1:3), paste0("S", 1:2))
    )
    sa <- data.frame(
        FullRunName = colnames(m_pep),
        Lab = c("A", "A"), Condition = c("X", "Y"),
        row.names = colnames(m_pep),
        stringsAsFactors = FALSE
    )

    pbf <- ProBatchFeatures(
        data_matrix       = m_pep,
        sample_annotation = sa,
        sample_id_col     = "FullRunName",
        level             = "peptide" # => "peptide::raw"
    )

    # Add a child assay with IDENTICAL rownames (1:1). Different pipeline.
    pbf2 <- pb_add_level(
        object      = pbf,
        from        = "peptide::raw",
        new_matrix  = m_pep,
        to_level    = "peptide",
        to_pipeline = "copy"
    )

    # The child exists
    expect_true("peptide::copy" %in% names(pbf2))

    # One-to-one link should be established (hits == nrow)
    al <- QFeatures::assayLink(pbf2, "peptide::copy")
    expect_s4_class(al, "AssayLink")
    expect_equal(length(al@hits), nrow(m_pep))

    # colData preserved and ordered
    se_from <- pbf2[["peptide::raw"]]
    se_to <- pbf2[["peptide::copy"]]
    expect_identical(
        rownames(SummarizedExperiment::colData(se_from)),
        rownames(SummarizedExperiment::colData(se_to))
    )

    # Log updated by pb_add_level
    log <- get_operation_log(pbf2)
    expect_true(nrow(log) >= 1L)
    expect_true(any(grepl("^add_level\\(peptide\\)_", as.character(log$step))))
})

test_that("pb_add_level: many-to-one linking (peptide -> protein) with map_strategy=first/longest (no splitting)", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    # Peptide matrix (3 peptides x 2 samples)
    m_pep <- matrix(
        c(
            10, 11,
            20, 21,
            30, 31
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(c("PEP1", "PEP2", "PEP3"), c("S1", "S2"))
    )
    sa <- data.frame(
        FullRunName = c("S1", "S2"),
        Lab = "LAB", Condition = c("C1", "C2"),
        row.names = c("S1", "S2"),
        stringsAsFactors = FALSE
    )
    pbf <- ProBatchFeatures(m_pep, sa, "FullRunName", level = "peptide")

    # Protein matrix (child rows are exact group strings)
    m_prot <- matrix(
        c(
            100, 101,
            200, 201,
            300, 301,
            400, 401
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(c("A", "B", "BBB;B", "C"), c("S1", "S2"))
    )

    # Mapping: duplicates per parent possible; selection among EXACT child strings only
    map <- data.frame(
        Precursor.Id = c("PEP1", "PEP2", "PEP3", "PEP3", "PEP4"),
        Protein.Ids = c("A", "A", "B", "BBB;B", "A"),
        stringsAsFactors = FALSE
    )

    # -------- Strategy: first (PEP3 -> "B" because it's the first exact candidate)
    p_first <- pb_add_level(
        object = pbf,
        from = "peptide::raw",
        new_matrix = m_prot,
        to_level = "protein",
        mapping_df = map,
        from_id = "Precursor.Id",
        to_id = "Protein.Ids",
        map_strategy = "first"
    )
    rf <- SummarizedExperiment::rowData(p_first[["peptide::raw"]])$ProteinID
    expect_identical(unname(rf), c("A", "A", "B"))

    al_first <- QFeatures::assayLink(p_first, "protein::raw")
    expect_s4_class(al_first, "AssayLink")
    expect_true(length(al_first@hits) >= 3L)

    # -------- Strategy: longest (PEP3 -> "BBB;B" : 2 proteins vs 1 in "B"; tie-breaks then by nchar, then order)
    p_longest <- pb_add_level(
        object = pbf,
        from = "peptide::raw",
        new_matrix = m_prot,
        to_level = "protein",
        mapping_df = map,
        from_id = "Precursor.Id",
        to_id = "Protein.Ids",
        map_strategy = "longest"
    )
    rl <- SummarizedExperiment::rowData(p_longest[["peptide::raw"]])$ProteinID
    expect_identical(unname(rl), c("A", "A", "BBB;B"))

    al_long <- QFeatures::assayLink(p_longest, "protein::raw")
    expect_s4_class(al_long, "AssayLink")
    expect_true(length(al_long@hits) >= 3L)

    # colData preserved in both
    expect_identical(
        rownames(SummarizedExperiment::colData(p_first[["peptide::raw"]])),
        rownames(SummarizedExperiment::colData(p_first[["protein::raw"]]))
    )
    expect_identical(
        rownames(SummarizedExperiment::colData(p_longest[["peptide::raw"]])),
        rownames(SummarizedExperiment::colData(p_longest[["protein::raw"]]))
    )
})

test_that("pb_add_level: map_strategy='as_is' errors on duplicate parents", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    m_pep <- matrix(1,
        nrow = 2, ncol = 2,
        dimnames = list(c("PEP1", "PEP2"), c("S1", "S2"))
    )
    sa <- data.frame(
        FullRunName = c("S1", "S2"),
        Lab = "L", Condition = c("C1", "C2"),
        row.names = c("S1", "S2"), stringsAsFactors = FALSE
    )
    pbf <- ProBatchFeatures(m_pep, sa, "FullRunName", level = "peptide")

    m_prot <- matrix(1,
        nrow = 2, ncol = 2,
        dimnames = list(c("A", "B"), c("S1", "S2"))
    )
    # Duplicate mapping for PEP1 -> A and PEP1 -> B
    map_dup <- data.frame(
        Precursor.Id = c("PEP1", "PEP1", "PEP2"),
        Protein.Ids = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )

    expect_error(
        pb_add_level(
            object = pbf,
            from = "peptide::raw",
            new_matrix = m_prot,
            to_level = "protein",
            mapping_df = map_dup,
            from_id = "Precursor.Id",
            to_id = "Protein.Ids",
            map_strategy = "as_is"
        ),
        "map_strategy='as_is'"
    )
})

test_that(".pb_add_assay_with_link requires identical row order for 1:1 linking (no set-equality)", {
    skip_if_not_installed("QFeatures")
    skip_if_not_installed("SummarizedExperiment")

    set.seed(42)
    m <- matrix(round(rnorm(12), 3),
        nrow = 4,
        dimnames = list(paste0("f", 1:4), paste0("S", 1:3)[1:3][1:3][1:3])
    )
    sa <- data.frame(
        FullRunName = colnames(m),
        Lab = "L", Condition = c("C1", "C2", "C3")[1:ncol(m)],
        row.names = colnames(m), stringsAsFactors = FALSE
    )
    pbf <- ProBatchFeatures(m, sa, "FullRunName", name = "raw")

    # Create a new SE with shuffled row order
    m2 <- m[sample(rownames(m)), ]
    se_new <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = m2),
        colData = SummarizedExperiment::colData(pbf[["feature::raw"]])
    )
    pbf2 <- proBatch:::.pb_add_assay_with_link(pbf,
        se = se_new,
        to = "feature::copy", from = "feature::raw"
    )

    # With setequal() requirement, shuffled order MUST link
    expect_true(isTRUE(S4Vectors::metadata(pbf2)$linked_last))
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
    se1 <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = m),
        colData = colData(pbf) # object-level colData (good)
    )
    pbf <- proBatch:::.pb_add_assay_with_link(pbf, se = se1, to = "feature::copy", from = "feature::raw")
    expect_true(isTRUE(S4Vectors::metadata(pbf)$linked_last))

    # Now create a NEW assay but (intentionally) give it assay-level colData copy and mutate DateTime
    se2 <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = m),
        colData = SummarizedExperiment::colData(pbf[["feature::copy"]])
    )
    SummarizedExperiment::colData(se2)$DateTime[1] <- se2$DateTime[1] + 3600 # conflict

    SummarizedExperiment::colData(se2) <- colData(pbf)
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
    se2 <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = m2),
        colData = colData(pbf)
    )
    pbf2 <- proBatch:::.pb_add_assay_with_link(pbf, se = se2, to = "feature::copy", from = "feature::raw")
    expect_true(isTRUE(S4Vectors::metadata(pbf2)$linked_last))
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
    expect_identical(rownames(colData(pbf2)), rownames(SummarizedExperiment::colData(pbf2[["protein::raw"]])))
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
