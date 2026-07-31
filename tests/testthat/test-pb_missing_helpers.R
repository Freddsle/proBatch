testthat::skip_if_not_installed("QFeatures")

make_test_pbf <- function(mat, sa_extra = NULL) {
    stopifnot(!is.null(colnames(mat)))
    sa <- data.frame(Sample = colnames(mat), stringsAsFactors = FALSE)
    if (!is.null(sa_extra)) {
        stopifnot(is.data.frame(sa_extra))
        if (nrow(sa_extra) != nrow(sa)) {
            stop("`sa_extra` must have one row per sample.")
        }
        sa <- cbind(sa, sa_extra)
    }
    suppressMessages(ProBatchFeatures(
        data_matrix = mat,
        sample_annotation = sa,
        sample_id_col = "Sample",
        name = "raw",
        level = "feature"
    ))
}

run_missing_filter_retry <- function(
  kind, object, source_name, target_name, params = list()
) {
    args <- list(
        object = object,
        pbf_name = source_name,
        inplace = FALSE,
        final_name = target_name
    )
    if (identical(kind, "group")) {
        args$group_cols <- "Batch"
    }
    do.call(
        if (identical(kind, "filter")) pb_filterNA else pb_groupfilterNA,
        c(args, params)
    )
}

virtualise_test_assay <- function(object, source_name) {
    suppressMessages(suppressWarnings(
        object[, , source_name, drop = FALSE]
    ))
}

test_that("pb_zeroIsNA converts zeros to missing values and logs the step", {
    mat <- matrix(
        c(0, 1, 2, 0, 3, 4),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)
    zero_count <- sum(assay(pbf[[assay_name]], "intensity") == 0)
    expect_gt(zero_count, 0L)

    res <- suppressMessages(pb_zeroIsNA(pbf))
    expect_s4_class(res, "ProBatchFeatures")

    updated <- assay(res[[assay_name]], "intensity")
    expect_equal(sum(is.na(updated)), zero_count)

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "zeroIsNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]], list())
})

test_that("pb_infIsNA replaces infinities and records the operation", {
    mat <- matrix(
        c(Inf, 2, 3, 4, 5, Inf),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)
    inf_count <- sum(is.infinite(assay(pbf[[assay_name]], "intensity")))
    expect_gt(inf_count, 0L)

    res <- suppressMessages(pb_infIsNA(pbf))
    expect_s4_class(res, "ProBatchFeatures")

    updated <- assay(res[[assay_name]], "intensity")
    expect_equal(sum(is.na(updated)), inf_count)

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "infIsNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]], list())
})

test_that("pb_nNA returns per-assay results for multiple assays", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, NA),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)

    # Create an additional stored assay with different missing pattern
    se_dup <- pbf[[assay_name]]
    alt_mat <- assay(se_dup, "intensity")
    alt_mat[1, 1] <- NA
    assay(se_dup, "intensity") <- alt_mat
    new_name <- paste0(assay_name, "_alt")
    pbf <- addAssay(pbf, se_dup, name = new_name)

    res <- pb_nNA(pbf, c(assay_name, new_name))
    expect_type(res, "list")
    expect_identical(names(res), c(assay_name, new_name, "nNA"))

    qf1 <- QFeatures(setNames(list(pbf[[assay_name]]), assay_name))
    qf2 <- QFeatures(setNames(list(pbf[[new_name]]), new_name))

    expect_type(res[[1]], "list")
    expect_type(res[[2]], "list")

    expect_identical(
        res[[1]],
        nNA(qf1, i = assay_name)
    )
    expect_identical(
        res[[2]],
        nNA(qf2, i = new_name)
    )
})

test_that("pb_filterNA stores filtered assays when not operating in place", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    original_names <- names(pbf)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_filterNA(pbf, inplace = FALSE))
    expect_s4_class(res, "ProBatchFeatures")

    expect_setequal(original_names, names(pbf))
    expect_equal(length(names(res)), length(original_names) + 1L)

    new_name <- setdiff(names(res), original_names)
    expect_length(new_name, 1L)
    expect_match(new_name, paste0(assay_name, "_filteredNA"))

    filtered <- assay(res[[new_name]], "intensity")
    expect_true(is.matrix(filtered))
    expect_false(anyNA(filtered))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "filterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), new_name)
    expect_true(isFALSE(log$params[[1]]$inplace) || is.null(log$params[[1]]$inplace))
})

test_that("pb_filterNA modifies stored assays in place when requested", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_filterNA(pbf, inplace = TRUE))
    expect_s4_class(res, "ProBatchFeatures")

    expect_identical(names(res), names(pbf))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_true(is.matrix(filtered))
    expect_false(anyNA(filtered))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "filterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]]$inplace, TRUE)
})

test_that("pb_filterNA preserves and prunes links during in-place filtering", {
    mat <- matrix(
        c(1, NA, 3, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    source_name <- pb_current_assay(pbf)
    linked_name <- "feature::copy"
    pbf <- proBatch:::.pb_add_assay_with_link(
        pbf,
        se = pbf[[source_name]],
        to = linked_name,
        from = source_name
    )

    result <- expect_no_warning(
        suppressMessages(
            pb_filterNA(
                pbf,
                pbf_name = linked_name,
                inplace = TRUE,
                pNA = 0
            )
        )
    )

    expect_true(validObject(result))
    expect_identical(rownames(result[[linked_name]]), c("f2", "f3"))
    link <- QFeatures::assayLink(result, linked_name)
    expect_identical(link@from, source_name)
    expect_identical(
        as.character(as.data.frame(link@hits)$names_to),
        c("f2", "f3")
    )
})

test_that("pb_filterNA validates final_name length when creating new assays", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)

    expect_error(
        pb_filterNA(pbf, inplace = FALSE, final_name = c("a", "b")),
        "`final_name` must contain exactly one name per selected assay"
    )
})

test_that("missing filters validate explicit final_name values consistently", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep("B1", ncol(mat)))
    )
    assay_name <- pb_current_assay(pbf)

    for (bad_name in list(NA_character_, "", 1)) {
        expected <- if (is.character(bad_name)) {
            "`final_name` must contain only non-missing, non-empty names."
        } else {
            "`final_name` must be a character vector."
        }

        expect_error(
            suppressMessages(pb_filterNA(
                pbf,
                pbf_name = assay_name,
                final_name = bad_name
            )),
            expected,
            fixed = TRUE
        )
        expect_error(
            suppressMessages(pb_groupfilterNA(
                pbf,
                pbf_name = assay_name,
                group_cols = "Batch",
                final_name = bad_name
            )),
            expected,
            fixed = TRUE
        )
    }
})

test_that("missing filters reject non-exact explicit stored-assay collisions", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep("B1", ncol(mat)))
    )
    assay_name <- pb_current_assay(pbf)

    expect_error(
        suppressMessages(pb_filterNA(
            pbf,
            pbf_name = assay_name,
            final_name = assay_name
        )),
        "different parent or stable lineage origin",
        fixed = TRUE
    )
    expect_error(
        suppressMessages(pb_groupfilterNA(
            pbf,
            pbf_name = assay_name,
            group_cols = "Batch",
            final_name = assay_name
        )),
        "different parent or stable lineage origin",
        fixed = TRUE
    )
    expect_identical(names(pbf), assay_name)
    expect_equal(nrow(get_operation_log(pbf)), 0L)
})

test_that("missing filters reject non-exact explicit virtual-target collisions", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep("B1", ncol(mat)))
    )
    assay_name <- pb_current_assay(pbf)
    virtual_name <- "reserved_missing_result"
    pbf <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "log2",
        fun = "log_transform_dm",
        from = assay_name,
        to = virtual_name,
        params = list(log_base = 2, offset = 1)
    )
    expect_false(virtual_name %in% names(pbf))

    expect_error(
        suppressMessages(pb_filterNA(
            pbf,
            pbf_name = assay_name,
            final_name = virtual_name
        )),
        "different parent or stable lineage origin",
        fixed = TRUE
    )
    expect_error(
        suppressMessages(pb_groupfilterNA(
            pbf,
            pbf_name = assay_name,
            group_cols = "Batch",
            final_name = virtual_name
        )),
        "different parent or stable lineage origin",
        fixed = TRUE
    )
    expect_identical(
        as.character(get_operation_log(pbf)$to),
        virtual_name
    )
})

test_that("missing filters make stored and virtual exact retries idempotent", {
    mat <- matrix(
        c(
            1, NA, 3, 4,
            5, 6, 7, 8
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:4))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep(c("B1", "B2"), each = 2))
    )
    source_name <- pb_current_assay(pbf)

    for (kind in c("filter", "group")) {
        target_name <- paste0(kind, "_retry")
        stored <- suppressMessages(run_missing_filter_retry(
            kind,
            pbf,
            source_name,
            target_name
        ))
        stored_data <- assay(stored[[target_name]], "intensity")
        stored_log <- get_operation_log(stored)

        stored_retry <- suppressMessages(run_missing_filter_retry(
            kind,
            stored,
            source_name,
            target_name
        ))
        expect_identical(names(stored_retry), names(stored))
        expect_identical(
            assay(stored_retry[[target_name]], "intensity"),
            stored_data
        )
        expect_identical(get_operation_log(stored_retry), stored_log)
        expect_identical(
            get_operation_log(stored_retry)$timestamp,
            stored_log$timestamp
        )

        virtual <- virtualise_test_assay(stored, source_name)
        expect_false(target_name %in% names(virtual))
        expect_identical(get_operation_log(virtual), stored_log)

        virtual_retry <- suppressMessages(run_missing_filter_retry(
            kind,
            virtual,
            source_name,
            target_name
        ))
        expect_true(target_name %in% names(virtual_retry))
        expect_identical(
            assay(virtual_retry[[target_name]], "intensity"),
            stored_data
        )
        expect_identical(get_operation_log(virtual_retry), stored_log)
        expect_identical(
            get_operation_log(virtual_retry)$timestamp,
            stored_log$timestamp
        )
    }
})

test_that("missing-filter retries compare data, parameters, and parents", {
    mat <- matrix(
        seq_len(8),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:4))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep(c("B1", "B2"), each = 2))
    )
    source_name <- pb_current_assay(pbf)
    alternate_name <- "feature::alternate"
    pbf <- addAssay(pbf, pbf[[source_name]], name = alternate_name)
    cases <- list(
        filter = list(first = list(pNA = 0), retry = list(pNA = 0.5)),
        group = list(
            first = list(min_valid = 1L),
            retry = list(min_valid = 2L)
        )
    )

    for (kind in names(cases)) {
        data_target <- paste0(kind, "_data_conflict")
        data_stored <- suppressMessages(run_missing_filter_retry(
            kind,
            pbf,
            source_name,
            data_target
        ))
        changed_source <- data_stored[[source_name]]
        changed_data <- assay(changed_source, "intensity")
        changed_data["f2", "s1"] <- changed_data["f2", "s1"] + 100
        assay(changed_source, "intensity") <- changed_data
        data_stored[[source_name]] <- changed_source

        expect_error(
            suppressMessages(run_missing_filter_retry(
                kind,
                data_stored,
                source_name,
                data_target
            )),
            "conflicting data",
            fixed = TRUE
        )

        parameter_target <- paste0(kind, "_parameter_conflict")
        parameter_stored <- suppressMessages(run_missing_filter_retry(
            kind,
            pbf,
            source_name,
            parameter_target,
            params = cases[[kind]]$first
        ))
        expect_identical(
            assay(parameter_stored[[parameter_target]], "intensity"),
            mat
        )
        expect_error(
            suppressMessages(run_missing_filter_retry(
                kind,
                parameter_stored,
                source_name,
                parameter_target,
                params = cases[[kind]]$retry
            )),
            "different parent or stable lineage origin",
            fixed = TRUE
        )

        parent_target <- paste0(kind, "_parent_conflict")
        parent_stored <- suppressMessages(run_missing_filter_retry(
            kind,
            pbf,
            source_name,
            parent_target
        ))
        expect_error(
            suppressMessages(run_missing_filter_retry(
                kind,
                parent_stored,
                alternate_name,
                parent_target
            )),
            "different parent or stable lineage origin",
            fixed = TRUE
        )
    }
})

test_that("missing filters require unique explicit names for multiple assays", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep("B1", ncol(mat)))
    )
    assay_name <- pb_current_assay(pbf)
    alt_name <- paste0(assay_name, "_alt")
    pbf <- addAssay(pbf, pbf[[assay_name]], name = alt_name)
    selected <- c(assay_name, alt_name)

    expect_error(
        suppressMessages(pb_filterNA(
            pbf,
            pbf_name = selected,
            final_name = "filtered"
        )),
        "scalar names are not recycled",
        fixed = TRUE
    )
    expect_error(
        suppressMessages(pb_groupfilterNA(
            pbf,
            pbf_name = selected,
            group_cols = "Batch",
            final_name = "group_filtered"
        )),
        "scalar names are not recycled",
        fixed = TRUE
    )

    for (duplicated_names in list(
        c("filtered", "filtered"),
        c("group_filtered", "group_filtered")
    )) {
        expect_error(
            suppressMessages(pb_filterNA(
                pbf,
                pbf_name = selected,
                final_name = duplicated_names
            )),
            "`final_name` entries must be unique",
            fixed = TRUE
        )
        expect_error(
            suppressMessages(pb_groupfilterNA(
                pbf,
                pbf_name = selected,
                group_cols = "Batch",
                final_name = duplicated_names
            )),
            "`final_name` entries must be unique",
            fixed = TRUE
        )
    }

    filter_names <- c("filtered_first", "filtered_second")
    filtered <- suppressMessages(pb_filterNA(
        pbf,
        pbf_name = selected,
        final_name = filter_names
    ))
    expect_true(all(filter_names %in% names(filtered)))
    expect_identical(
        tail(as.character(get_operation_log(filtered)$to), 2L),
        filter_names
    )

    group_names <- c("group_filtered_first", "group_filtered_second")
    group_filtered <- suppressMessages(pb_groupfilterNA(
        pbf,
        pbf_name = selected,
        group_cols = "Batch",
        final_name = group_names
    ))
    expect_true(all(group_names %in% names(group_filtered)))
    expect_identical(
        tail(as.character(get_operation_log(group_filtered)$to), 2L),
        group_names
    )
})

test_that("missing filters disambiguate generated names across the result namespace", {
    mat <- matrix(
        c(1, NA, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = rep("B1", ncol(mat)))
    )
    assay_name <- pb_current_assay(pbf)
    stored_default <- paste0(assay_name, "_filteredNA")
    virtual_default <- paste0(assay_name, "_groupfilteredNA")
    pbf <- addAssay(pbf, pbf[[assay_name]], name = stored_default)
    pbf <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "log2",
        fun = "log_transform_dm",
        from = assay_name,
        to = virtual_default,
        params = list(log_base = 2, offset = 1)
    )

    filtered <- suppressMessages(pb_filterNA(
        pbf,
        pbf_name = assay_name
    ))
    expect_true(paste0(stored_default, ".1") %in% names(filtered))

    filtered_again <- suppressMessages(pb_filterNA(
        filtered,
        pbf_name = assay_name
    ))
    expect_true(paste0(stored_default, ".2") %in% names(filtered_again))

    group_filtered <- suppressMessages(pb_groupfilterNA(
        pbf,
        pbf_name = assay_name,
        group_cols = "Batch"
    ))
    expect_true(paste0(virtual_default, ".1") %in% names(group_filtered))
})

test_that("pb_groupfilterNA retains union of group-wise valid features in place", {
    mat <- matrix(
        c(
            1, 2, 3, 4,
            NA, NA, 5, 6,
            1, NA, 2, NA
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa_extra <- data.frame(
        Batch = rep(c("B1", "B2"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        inplace = TRUE
    ))

    expect_s4_class(res, "ProBatchFeatures")
    expect_identical(names(res), names(pbf))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_true(is.matrix(filtered))
    expect_identical(rownames(filtered), c("f1", "f2"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "groupfilterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), assay_name)
    expect_identical(log$params[[1]]$group_cols, "Batch")
    expect_identical(log$params[[1]]$min_valid, 2L)
    expect_null(log$params[[1]]$pNA)
})

test_that("pb_groupfilterNA reports links removed by in-place replacement", {
    mat <- matrix(
        c(1, 2, NA, 4, 5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(
        mat,
        sa_extra = data.frame(Batch = c("B1", "B1"))
    )
    source_name <- pb_current_assay(pbf)
    linked_name <- "feature::copy"
    pbf <- proBatch:::.pb_add_assay_with_link(
        pbf,
        se = pbf[[source_name]],
        to = linked_name,
        from = source_name
    )

    result <- NULL
    expect_warning(
        result <- suppressMessages(
            pb_groupfilterNA(
                pbf,
                pbf_name = linked_name,
                group_cols = "Batch",
                min_valid = 2L,
                inplace = TRUE
            )
        ),
        "Links between assays were lost/removed during replacement",
        fixed = TRUE
    )

    expect_true(validObject(result))
    expect_identical(rownames(result[[linked_name]]), c("f1", "f3"))
    link <- QFeatures::assayLink(result, linked_name)
    expect_true(is.na(link@from))
    expect_length(link@hits, 0L)
})

test_that("pb_groupfilterNA stores union of group-wise valid features when not operating in place", {
    mat <- matrix(
        c(
            1, 2, 3, 4,
            NA, NA, 5, 6,
            1, NA, 2, NA
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa_extra <- data.frame(
        Batch = rep(c("B1", "B2"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)
    original_names <- names(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        inplace = FALSE,
        final_name = "filtered_group"
    ))

    expect_s4_class(res, "ProBatchFeatures")
    expect_setequal(original_names, names(pbf))
    expect_equal(length(names(res)), length(original_names) + 1L)

    new_name <- setdiff(names(res), original_names)
    expect_identical(new_name, "filtered_group")

    filtered <- assay(res[[new_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f2"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_identical(as.character(log$step), "groupfilterNA")
    expect_identical(as.character(log$from), assay_name)
    expect_identical(as.character(log$to), new_name)
    expect_identical(log$params[[1]]$group_cols, "Batch")
    expect_identical(log$params[[1]]$min_valid, 2L)
    expect_null(log$params[[1]]$pNA)
})

test_that("pb_groupfilterNA respects pNA thresholds when supplied", {
    mat <- matrix(
        c(
            1, 2, 3, 4,
            NA, 2, NA, 4,
            NA, NA, 5, 6
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
    )
    sa_extra <- data.frame(
        Batch = rep(c("B1", "B2"), each = 2),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = NULL,
        pNA = 0.25,
        inplace = TRUE
    ))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f3"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_null(log$params[[1]]$min_valid)
    expect_equal(log$params[[1]]$pNA, 0.25)
})

test_that("pb_groupfilterNA combines min_valid and pNA per group", {
    mat <- matrix(
        c(
            1, 2, 3, 4, 5, 6,
            NA, 2, 3, 7, 8, 9,
            NA, NA, 5, 6, NA, NA
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:3), paste0("s", 1:6))
    )
    sa_extra <- data.frame(
        Batch = c(rep("B1", 3), rep("B2", 3)),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        pNA = 0.25,
        inplace = TRUE
    ))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f2"))

    log <- get_operation_log(res)
    expect_equal(nrow(log), 1L)
    expect_equal(log$params[[1]]$min_valid, 2L)
    expect_equal(log$params[[1]]$pNA, 0.25)
})

test_that("pb_groupfilterNA validates presence of grouping columns", {
    mat <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)

    expect_error(
        pb_groupfilterNA(pbf, group_cols = "Batch"),
        "missing group column(s)",
        fixed = TRUE
    )
})

test_that("pb_groupfilterNA errors when a group has too few samples", {
    mat <- matrix(
        c(
            1, 2, 3,
            4, 5, NA
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:3))
    )
    sa_extra <- data.frame(
        Batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )
    pbf <- make_test_pbf(mat, sa_extra = sa_extra)

    expect_error(
        pb_groupfilterNA(pbf, group_cols = "Batch", min_valid = 2L),
        "requires at least 2",
        fixed = TRUE
    )
})

test_that("pb_missing helpers error on non-materialised assays with guidance", {
    mat <- matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(paste0("f", 1:2), paste0("s", 1:2))
    )
    pbf <- make_test_pbf(mat)
    assay_name <- pb_current_assay(pbf)

    expect_error(
        pb_zeroIsNA(pbf, pbf_name = "feature::not_there"),
        "are not stored in the object",
        fixed = TRUE
    )

    logged_only <- proBatch:::.pb_add_log_entry(
        pbf,
        step = "log2",
        fun = "log_transform_dm",
        from = assay_name,
        to = paste0(assay_name, "_log2"),
        params = list()
    )

    expect_error(
        pb_zeroIsNA(logged_only, pbf_name = paste0(assay_name, "_log2")),
        "store_fast_steps = TRUE",
        fixed = TRUE
    )
})


# ------------------------------------------------------------------
# pb_groupfilterNA: mask_failing parameter
# ------------------------------------------------------------------

.mask_mat <- matrix(
    c(
        1, 2, 3, 4,
        NA, 7, 5, 6,
        8, 9, NA, 10
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("f", 1:3), paste0("s", 1:4))
)
.mask_sa <- data.frame(
    Batch = rep(c("B1", "B2"), each = 2),
    stringsAsFactors = FALSE
)

test_that("pb_groupfilterNA masks values in failing groups by default", {
    pbf <- make_test_pbf(.mask_mat, sa_extra = .mask_sa)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        inplace = TRUE
    ))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f2", "f3"))

    expect_true(is.na(filtered["f2", "s1"]))
    expect_true(is.na(filtered["f2", "s2"]))
    expect_equal(filtered["f2", "s3"], 5)
    expect_equal(filtered["f2", "s4"], 6)

    expect_equal(filtered["f3", "s1"], 8)
    expect_equal(filtered["f3", "s2"], 9)
    expect_true(is.na(filtered["f3", "s3"]))
    expect_true(is.na(filtered["f3", "s4"]))

    expect_equal(as.numeric(filtered["f1", ]), c(1, 2, 3, 4))

    log <- get_operation_log(res)
    expect_identical(log$params[[1]]$mask_failing, TRUE)
})

test_that("pb_groupfilterNA keeps failing-group values when masking is disabled", {
    pbf <- make_test_pbf(.mask_mat, sa_extra = .mask_sa)
    assay_name <- pb_current_assay(pbf)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        group_cols = "Batch",
        min_valid = 2L,
        mask_failing = FALSE,
        inplace = TRUE
    ))

    filtered <- assay(res[[assay_name]], "intensity")
    expect_identical(rownames(filtered), c("f1", "f2", "f3"))

    expect_true(is.na(filtered["f2", "s1"]))
    expect_equal(filtered["f2", "s2"], 7)
    expect_equal(filtered["f2", "s3"], 5)
    expect_equal(filtered["f2", "s4"], 6)

    expect_equal(filtered["f3", "s1"], 8)
    expect_equal(filtered["f3", "s2"], 9)
    expect_true(is.na(filtered["f3", "s3"]))
    expect_equal(filtered["f3", "s4"], 10)

    log <- get_operation_log(res)
    expect_identical(log$params[[1]]$mask_failing, FALSE)
})

test_that("pb_groupfilterNA preserves positional inplace compatibility", {
    pbf <- make_test_pbf(.mask_mat, sa_extra = .mask_sa)

    res <- suppressMessages(pb_groupfilterNA(
        pbf,
        NULL,
        "Batch",
        2L,
        NULL,
        TRUE
    ))

    expect_identical(names(res), names(pbf))
    log <- get_operation_log(res)
    expect_identical(log$params[[1]]$inplace, TRUE)
    expect_identical(log$params[[1]]$mask_failing, TRUE)
})
