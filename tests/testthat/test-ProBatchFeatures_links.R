test_that("internal add-assay-with-link links one-to-one when rows match", {
    pb_test_load_example_data()

    pbf <- ProBatchFeatures(example_proteome_matrix, example_sample_annotation, "FullRunName", name = "raw")

    # create a new assay with identical rownames -> should link 1:1
    se_new <- SummarizedExperiment(
        assays  = list(intensity = pb_as_wide(pbf)),
        colData = colData(pbf[["feature::raw"]])
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
    colData(se_mod)$DateTime <- NULL

    pbf3 <- proBatch:::.pb_add_assay_with_link(pbf2, se = se_mod, to = "feature::raw_changed", from = "feature::raw_copy")
    expect_false(isTRUE(S4Vectors::metadata(pbf3)$linked_last))
})

test_that("show() prints chain and step count", {
    pb_test_load_example_data()

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
    expect_match(out1, "Processing chain:\\n +\\[1\\] log")
    expect_match(out1, "Steps logged: 1")
})

test_that("validity checks fail on malformed oplog", {
    pb_test_load_example_data()

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
        rownames(colData(se_from)),
        rownames(colData(se_to))
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
    rf <- rowData(p_first[["peptide::raw"]])$ProteinID
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
    rl <- rowData(p_longest[["peptide::raw"]])$ProteinID
    expect_identical(unname(rl), c("A", "A", "BBB;B"))

    al_long <- QFeatures::assayLink(p_longest, "protein::raw")
    expect_s4_class(al_long, "AssayLink")
    expect_true(length(al_long@hits) >= 3L)

    # colData preserved in both
    expect_identical(
        rownames(colData(p_first[["peptide::raw"]])),
        rownames(colData(p_first[["protein::raw"]]))
    )
    expect_identical(
        rownames(colData(p_longest[["peptide::raw"]])),
        rownames(colData(p_longest[["protein::raw"]]))
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
    se_new <- SummarizedExperiment(
        assays  = list(intensity = m2),
        colData = colData(pbf[["feature::raw"]])
    )
    pbf2 <- proBatch:::.pb_add_assay_with_link(pbf,
        se = se_new,
        to = "feature::copy", from = "feature::raw"
    )

    # With setequal() requirement, shuffled order MUST link
    expect_true(isTRUE(S4Vectors::metadata(pbf2)$linked_last))
})
