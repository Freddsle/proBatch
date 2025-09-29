test_that("find_duplicated_columns detects duplicates in data frames", {
    df <- data.frame(
        a = c(1, NA, 3, NA),
        b = c(1, NA, 3, NA),
        c = c(0, 2, 4, 6),
        d = c("x", "y", "z", "w"),
        e = c("x", "y", "z", "w"),
        stringsAsFactors = FALSE
    )

    duplicates <- find_duplicated_columns(df)
    pairs <- vapply(duplicates, function(x) paste(sort(x), collapse = ","), character(1))

    expect_setequal(pairs, c("a,b", "d,e"))
    expect_identical(find_duplicated_columns(df["a"]), list())
})

test_that("find_duplicated_columns assigns default names when absent", {
    unnamed <- data.frame(c(1, 2), c(1, 2))
    colnames(unnamed) <- NULL

    duplicates <- find_duplicated_columns(unnamed)
    expect_identical(duplicates, list(c("V1", "V2")))
})

test_that("find_duplicated_columns.ProBatchFeatures inspects sample metadata by default", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")

    example_sample_annotation$Diet_copy <- example_sample_annotation$Diet

    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    dup_pairs <- find_duplicated_columns(pbf)
    pair_strings <- vapply(dup_pairs, function(x) paste(sort(x), collapse = ","), character(1))
    expect_true("Diet,Diet_copy" %in% pair_strings)
})

test_that("find_duplicated_columns.ProBatchFeatures can inspect rowData and custom data frames", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    example_sample_annotation$Diet_copy <- example_sample_annotation$Diet

    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    assay_name <- pb_current_assay(pbf)
    peptide_anno <- example_peptide_annotation
    peptide_anno$Protein_copy <- peptide_anno$ProteinName
    idx <- match(rownames(example_proteome_matrix), peptide_anno$peptide_group_label)
    peptide_anno <- peptide_anno[idx, , drop = FALSE]
    rownames(peptide_anno) <- peptide_anno$peptide_group_label
    peptide_anno$peptide_group_label <- NULL
    SummarizedExperiment::rowData(pbf[[assay_name]]) <- S4Vectors::DataFrame(peptide_anno)

    row_pairs <- find_duplicated_columns(pbf, component = "rowData", assay = assay_name)
    row_strings <- vapply(row_pairs, function(x) paste(sort(x), collapse = ","), character(1))
    expect_true("ProteinName,Protein_copy" %in% row_strings)

    custom_df <- data.frame(x = 1:3, y = 1:3)
    custom_pairs <- find_duplicated_columns(pbf, df = custom_df)
    expect_identical(custom_pairs, list(c("x", "y")))
})

test_that("metadata_column_summary handles data frames", {
    df <- data.frame(
        a = c(1, 1, NA, 2),
        b = c(NA, NA, NA, NA),
        c = c("x", "y", "z", "z"),
        stringsAsFactors = FALSE
    )

    summary_sorted <- metadata_column_summary(df)
    expect_equal(summary_sorted$colname, c("b", "a", "c"))
    expect_equal(summary_sorted$n_unique[summary_sorted$colname == "b"], 0)
    expect_equal(summary_sorted$n_NA[summary_sorted$colname == "a"], 1)
    expect_equal(summary_sorted$pct_NA[summary_sorted$colname == "a"], 25)

    summary_unsorted <- metadata_column_summary(df, sort = FALSE)
    expect_equal(summary_unsorted$colname, c("a", "b", "c"))
})

test_that("metadata_column_summary.ProBatchFeatures summarises object metadata", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    example_sample_annotation$Diet_copy <- example_sample_annotation$Diet

    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    col_summary <- metadata_column_summary(pbf)
    expect_true("Diet_copy" %in% col_summary$colname)
    diet_row <- col_summary[col_summary$colname == "Diet_copy", , drop = FALSE]
    expect_equal(diet_row$n_unique, length(unique(example_sample_annotation$Diet)))

    assay_name <- pb_current_assay(pbf)
    peptide_anno <- example_peptide_annotation
    peptide_anno$Protein_copy <- peptide_anno$ProteinName
    idx <- match(rownames(example_proteome_matrix), peptide_anno$peptide_group_label)
    peptide_anno <- peptide_anno[idx, , drop = FALSE]
    rownames(peptide_anno) <- peptide_anno$peptide_group_label
    peptide_anno$peptide_group_label <- NULL
    SummarizedExperiment::rowData(pbf[[assay_name]]) <- S4Vectors::DataFrame(peptide_anno)

    row_summary <- metadata_column_summary(pbf, component = "rowData", assay = assay_name)
    expect_true("Protein_copy" %in% row_summary$colname)

    custom_df <- data.frame(u = c(1, 2, 3), v = c(1, 2, 3))
    custom_summary <- metadata_column_summary(pbf, df = custom_df)
    expect_identical(custom_summary$colname, c("u", "v"))
})

test_that("filter_metadata_columns handles duplicates by pattern", {
    df <- data.frame(
        check.names = FALSE,
        `factor value intensity` = c(1, 2, 3),
        intensity = c(1, 2, 3),
        other = c(3, 2, 1)
    )

    filtered <- filter_metadata_columns(
        df,
        duplicate_keep = "pattern",
        duplicate_pattern = "^factor value"
    )

    expect_equal(colnames(filtered), c("factor value intensity", "other"))
    dropped <- attr(filtered, "dropped_duplicates")
    expect_identical(dropped, "intensity")

    filtered_last <- filter_metadata_columns(df, duplicate_keep = "last")
    expect_equal(colnames(filtered_last), c("intensity", "other"))
})

test_that("filter_metadata_columns removes columns based on missingness thresholds", {
    df <- data.frame(
        a = c(1, NA, 3, 4),
        b = c(NA, NA, NA, NA),
        c = c(1, 2, NA, NA),
        stringsAsFactors = FALSE
    )

    filtered <- filter_metadata_columns(df, min_non_na = 3)
    expect_equal(colnames(filtered), "a")
    expect_equal(attr(filtered, "dropped_missing"), c("b", "c"))

    filtered_pct <- filter_metadata_columns(df, max_pct_na = 50)
    expect_equal(colnames(filtered_pct), "a")
})

test_that("filter_metadata_columns.ProBatchFeatures filters metadata and can update inplace", {
    data(example_proteome_matrix, package = "proBatch")
    data(example_sample_annotation, package = "proBatch")
    data(example_peptide_annotation, package = "proBatch")

    example_sample_annotation$Diet_copy <- example_sample_annotation$Diet

    pbf <- ProBatchFeatures(
        data_matrix = example_proteome_matrix,
        sample_annotation = example_sample_annotation,
        sample_id_col = "FullRunName",
        name = "raw"
    )

    filtered_df <- filter_metadata_columns(
        pbf,
        duplicate_keep = "pattern",
        duplicate_pattern = "_copy$"
    )
    expect_false("Diet" %in% colnames(filtered_df))
    expect_true("Diet_copy" %in% colnames(filtered_df))

    assay_name <- pb_current_assay(pbf)
    peptide_anno <- example_peptide_annotation
    peptide_anno$Protein_copy <- peptide_anno$ProteinName
    idx <- match(rownames(example_proteome_matrix), peptide_anno$peptide_group_label)
    peptide_anno <- peptide_anno[idx, , drop = FALSE]
    rownames(peptide_anno) <- peptide_anno$peptide_group_label
    peptide_anno$peptide_group_label <- NULL
    SummarizedExperiment::rowData(pbf[[assay_name]]) <- S4Vectors::DataFrame(peptide_anno)

    updated <- filter_metadata_columns(
        pbf,
        component = "rowData",
        duplicate_keep = "pattern",
        duplicate_pattern = "^Protein_copy$",
        inplace = TRUE
    )

    row_meta <- SummarizedExperiment::rowData(updated[[assay_name]])
    expect_false("ProteinName" %in% colnames(row_meta))
    expect_true("Protein_copy" %in% colnames(row_meta))
    info <- metadata(updated)$filter_metadata_columns
    expect_type(info, "list")
    expect_identical(info$component, "rowData")
    expect_identical(info$assay, assay_name)
    expect_true("ProteinName" %in% info$dropped_columns)
})
