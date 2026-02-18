mock_intragroup_plot <- function(tag) { # helper returning simple ggplot with stored data
    df <- data.frame(
        Normalization = tag,
        value = seq_along(tag)
    )
    ggplot2::ggplot(df, ggplot2::aes(x = Normalization, y = value)) +
        ggplot2::geom_col()
}

test_that("plot_intragroup_variation.default forwards matrices to PRONE", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(10, 20, 30, 40, 50, 60),
        nrow = 2,
        dimnames = list(
            c("feat1", "feat2"),
            c("s1", "s2", "s3")
        )
    )
    sample_ann <- data.frame(
        FullRunName = c("s3", "s1", "s2"),
        Condition = c("A", "A", "B"),
        stringsAsFactors = FALSE
    )

    tmp_dir <- file.path(tempdir(), "intragr_default")
    if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$cor <- list(se = se, ain = ain, condition = condition, method = method)
            mock_intragroup_plot("cor")
        },
        plot_intragroup_PCV = function(se, ain, condition, diff) {
            captured$pcv <- list(se = se, ain = ain, condition = condition, diff = diff)
            mock_intragroup_plot("pcv")
        },
        .package = "PRONE"
    )

    res <- plot_intragroup_variation(
        dm,
        sample_annotation = sample_ann,
        group_col = "Condition",
        metrics = "correlation",
        correlation_method = "spearman",
        plot_title = "peptide::raw",
        path_to_save_results = tmp_dir
    )

    expect_s3_class(res, "ggplot")
    metrics_attr <- attr(res, "pb_intragroup_metrics")
    expect_true("correlation" %in% names(metrics_attr))
    saved_attr <- attr(res, "pb_intragroup_saved_files")
    expect_true("correlation" %in% names(saved_attr))
    expect_true(all(file.exists(unlist(saved_attr))))

    expect_s4_class(captured$cor$se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::assayNames(captured$cor$se), "peptide::raw")
    expect_equal(captured$cor$ain, "peptide::raw")
    expect_equal(captured$cor$condition, "Condition")
    expect_equal(SummarizedExperiment::colnames(captured$cor$se), colnames(dm))
    expect_true(isTRUE(captured$pcv$diff))
})

test_that("plot_intragroup_variation.default honors fill_the_missing", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(NA, 25, 50, 75, 100, 125),
        nrow = 2,
        dimnames = list(
            c("feat1", "feat2"),
            c("s1", "s2", "s3")
        )
    )
    sample_ann <- data.frame(
        FullRunName = c("s1", "s2", "s3"),
        Condition = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$assay <- SummarizedExperiment::assay(se, ain)
            mock_intragroup_plot("cor")
        },
        .package = "PRONE"
    )

    expect_warning(
        expect_warning(
            plot_intragroup_variation(
                dm,
                sample_annotation = sample_ann,
                group_col = "Condition",
                metrics = "correlation",
                fill_the_missing = 0
            ),
            "Intragroup diagnostics cannot operate with missing values",
            fixed = TRUE
        ),
        "filling missing values with 0",
        fixed = TRUE
    )

    expect_false(anyNA(captured$assay))
    expect_identical(captured$assay[1, 1], 0)
})

test_that("plot_intragroup_variation.default handles multiple grouping columns", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(5, 15, 25, 35, 45, 55),
        nrow = 2,
        dimnames = list(
            c("feat1", "feat2"),
            c("s1", "s2", "s3")
        )
    )
    sample_ann <- data.frame(
        FullRunName = c("s1", "s2", "s3"),
        Condition = c("A", "B", "A"),
        Instrument = c("I1", "I1", "I2"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$condition <- c(captured$condition, condition)
            mock_intragroup_plot(condition)
        },
        .package = "PRONE"
    )

    gg <- plot_intragroup_variation(
        dm,
        sample_annotation = sample_ann,
        group_col = c("Condition", "Instrument", "Condition"),
        metrics = "correlation",
        plot_title = "assay_one"
    )

    expect_s3_class(gg, "ggplot")
    expect_setequal(captured$condition, c("Condition", "Instrument"))
    expect_equal(sort(unique(as.character(gg$data$group_column))), sort(c("Condition", "Instrument")))
    expect_equal(gg$labels$x, "Grouping factor")
    expect_match(gg$labels$title, "^Intragroup correlation", perl = TRUE)
    expect_match(gg$labels$title, "\\nAssay: ")
})

test_that("plot_intragroup_variation.ProBatchFeatures iterates across assays", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        seq_len(8),
        nrow = 2,
        dimnames = list(
            c("feat1", "feat2"),
            paste0("sample", 1:4)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("sample", 1:4),
        Condition = rep(c("A", "B"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    )

    out_dir <- file.path(tempdir(), "intragr_pbf")
    if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$ain <- c(captured$ain, ain)
            mock_intragroup_plot("cor")
        },
        .package = "PRONE"
    )

    res <- plot_intragroup_variation(
        pbf,
        group_col = "Condition",
        metrics = "correlation",
        path_to_save_results = out_dir
    )

    expect_s3_class(res, "ggplot")
    expect_true(dir.exists(file.path(out_dir, "feature::raw")))
    saved_attr <- attr(res, "pb_intragroup_saved_files")
    expect_true("feature::raw" %in% names(saved_attr))
    expect_equal(unique(captured$ain), "feature::raw")
})

test_that("plot_intragroup_variation.ProBatchFeatures supports multiple metrics", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        seq_len(8),
        nrow = 2,
        dimnames = list(
            c("feat1", "feat2"),
            paste0("sample", 1:4)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("sample", 1:4),
        Condition = rep(c("A", "B"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    )

    captured <- new.env(parent = emptyenv())
    captured$calls <- character()
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$calls <- c(captured$calls, sprintf("correlation:%s", condition))
            mock_intragroup_plot("cor_multi")
        },
        plot_intragroup_PCV = function(se, ain, condition, diff) {
            captured$calls <- c(captured$calls, sprintf("PCV:%s", condition))
            mock_intragroup_plot("pcv_multi")
        },
        .package = "PRONE"
    )

    res <- plot_intragroup_variation(
        pbf,
        group_col = "Condition",
        metrics = c("correlation", "PCV")
    )

    expect_true(inherits(res, "ggplot") || grid::is.grob(res))
    expect_setequal(captured$calls, c("correlation:Condition", "PCV:Condition"))
    metrics_attr <- attr(res, "pb_intragroup_metrics")
    expect_true("feature::raw" %in% names(metrics_attr))
    expect_setequal(names(metrics_attr[["feature::raw"]]), c("correlation", "PCV"))
})

test_that("plot_intragroup_variation.ProBatchFeatures facets per grouping column", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        seq_len(12),
        nrow = 3,
        dimnames = list(
            paste0("feat", 1:3),
            paste0("sample", 1:4)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("sample", 1:4),
        Condition = rep(c("A", "B"), each = 2),
        Instrument = rep(c("MS1", "MS2"), times = 2),
        stringsAsFactors = FALSE
    )
    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    pbf <- suppressMessages(pb_transform(
        pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    ))

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$condition <- c(captured$condition, condition)
            captured$ain <- c(captured$ain, ain)
            mock_intragroup_plot(condition)
        },
        .package = "PRONE"
    )

    gg <- plot_intragroup_variation(
        pbf,
        group_col = c("Condition", "Instrument"),
        metrics = "correlation",
        plot_ncol = 1
    )

    expect_s3_class(gg, "ggplot")
    expect_equal(sort(unique(as.character(gg$data$group_column))), sort(c("Condition", "Instrument")))
    expect_equal(length(unique(as.character(gg$data$assay_display))), length(names(pbf)))
    expect_equal(sort(unique(captured$condition)), sort(c("Condition", "Instrument")))
    expect_equal(sort(unique(captured$ain)), sort(names(pbf)))
    expect_s3_class(gg$facet, "FacetWrap")
    expect_equal(gg$facet$params$ncol, 1)
    facet_var <- rlang::as_label(gg$facet$params$facets[[1]])
    expect_equal(facet_var, "group_column")
    expect_equal(gg$labels$x, "Assay")
    expect_equal(gg$labels$title, "Intragroup correlation")
})

test_that("plot_intragroup_variation.ProBatchFeatures honors explicit pbf_name", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        seq_len(12),
        nrow = 3,
        dimnames = list(
            paste0("feat", 1:3),
            paste0("sample", 1:4)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("sample", 1:4),
        Condition = rep(c("A", "B"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- suppressMessages(ProBatchFeatures(
        data_matrix = dm,
        sample_annotation = sample_ann,
        sample_id_col = "FullRunName",
        name = "feature::raw"
    ))
    pbf <- suppressMessages(pb_transform(
        pbf,
        from = "feature::raw",
        steps = "log2",
        store_fast_steps = TRUE
    ))

    selected_assay <- "feature::raw"

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$ain <- c(captured$ain, ain)
            mock_intragroup_plot("cor")
        },
        .package = "PRONE"
    )

    gg <- plot_intragroup_variation(
        pbf,
        pbf_name = selected_assay,
        group_col = "Condition",
        metrics = "correlation"
    )

    expect_s3_class(gg, "ggplot")
    expect_equal(unique(captured$ain), selected_assay)
    expect_equal(unique(as.character(gg$data$assay_display)), selected_assay)
})

test_that("plot_intragroup_variation.default validates intragroup column", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(1:6),
        nrow = 2,
        dimnames = list(
            paste0("f", 1:2),
            paste0("s", 1:3)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("s", c(3, 1, 2)),
        Condition = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )

    expect_error(
        plot_intragroup_variation(
            dm,
            sample_annotation = sample_ann,
            group_col = "MissingCol"
        ),
        "group_col value\\(s\\) not found",
        fixed = FALSE
    )

    sample_ann$Condition[1] <- NA
    expect_error(
        plot_intragroup_variation(
            dm,
            sample_annotation = sample_ann,
            group_col = "Condition"
        ),
        "contains missing values",
        fixed = TRUE
    )
})

test_that("plot_intragroup_variation saves sanitized metric tables", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        seq_len(8),
        nrow = 2,
        dimnames = list(
            c("feat1", "feat2"),
            paste0("s", 1:4)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("s", c(2, 1, 4, 3)),
        Condition = rep(c("A", "B"), each = 2),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$last_ain <- ain
            mock_intragroup_plot("cor_save")
        },
        .package = "PRONE"
    )

    save_dir <- file.path(withr::local_tempdir(), "metrics output")
    res <- plot_intragroup_variation(
        dm,
        sample_annotation = sample_ann,
        group_col = "Condition",
        metrics = "correlation",
        path_to_save_results = save_dir,
        plot_title = "peptide::raw/Normalized Step"
    )

    expect_true(dir.exists(save_dir))
    saved_attr <- attr(res, "pb_intragroup_saved_files")
    expect_named(saved_attr, "correlation")
    saved_file <- saved_attr$correlation
    expect_true(file.exists(saved_file))
    expect_equal(
        basename(saved_file),
        "intragroup_correlation_peptide_raw_Normalized_Step.csv"
    )
    expect_equal(captured$last_ain, "peptide::raw/Normalized Step")
})

test_that("plot_intragroup_variation supports multiple metrics", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(1:6),
        nrow = 2,
        dimnames = list(
            paste0("f", 1:2),
            paste0("s", 1:3)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("s", c(3, 1, 2)),
        Condition = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    captured$calls <- character()
    testthat::local_mocked_bindings(
        plot_intragroup_correlation = function(se, ain, condition, method) {
            captured$calls <- c(captured$calls, sprintf("correlation:%s", condition))
            mock_intragroup_plot("cor_multi")
        },
        plot_intragroup_PCV = function(se, ain, condition, diff) {
            captured$calls <- c(captured$calls, sprintf("PCV:%s", condition))
            mock_intragroup_plot("pcv_multi")
        },
        .package = "PRONE"
    )

    res <- plot_intragroup_variation(
        dm,
        sample_annotation = sample_ann,
        group_col = "Condition",
        metrics = c("correlation", "PCV")
    )

    expect_true(inherits(res, "ggplot") || grid::is.grob(res))
    metrics_attr <- attr(res, "pb_intragroup_metrics")
    expect_setequal(names(metrics_attr), c("correlation", "PCV"))
    expect_setequal(captured$calls, c("correlation:Condition", "PCV:Condition"))
    expect_null(attr(res, "pb_intragroup_saved_files"))
})

test_that("plot_intragroup_variation retries diff metrics with diff = FALSE when needed", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(1:6),
        nrow = 2,
        dimnames = list(
            paste0("f", 1:2),
            paste0("s", 1:3)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("s", c(3, 1, 2)),
        Condition = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    captured$diff <- logical()
    testthat::local_mocked_bindings(
        plot_intragroup_PCV = function(se, ain, condition, diff) {
            captured$diff <- c(captured$diff, diff)
            if (isTRUE(diff)) {
                stop("Log2 data not in SummarizedExperiment! Difference not applicable.")
            }
            mock_intragroup_plot("pcv_retry")
        },
        .package = "PRONE"
    )

    res <- NULL
    expect_warning(
        res <- plot_intragroup_variation(
            dm,
            sample_annotation = sample_ann,
            group_col = "Condition",
            metrics = "PCV",
            pcv_diff = TRUE
        ),
        "retrying with diff = FALSE",
        fixed = TRUE
    )

    expect_true(inherits(res, "ggplot") || grid::is.grob(res))
    expect_identical(captured$diff, c(TRUE, FALSE))
})

test_that("plot_intragroup_variation requires non-empty grouping columns", {
    testthat::skip_if_not_installed("PRONE")

    dm <- matrix(
        c(1:6),
        nrow = 2,
        dimnames = list(
            paste0("f", 1:2),
            paste0("s", 1:3)
        )
    )
    sample_ann <- data.frame(
        FullRunName = paste0("s", c(3, 1, 2)),
        Condition = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )

    expect_error(
        plot_intragroup_variation(
            dm,
            sample_annotation = sample_ann,
            group_col = ""
        ),
        "`group_col` must contain at least one non-empty column name.",
        fixed = TRUE
    )
})
