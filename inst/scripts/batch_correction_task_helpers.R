# Helper utilities for running one or many batch-correction tasks.
# Reused and adapted from inst/scripts/03_data_processing.R task workflow.

`%||%` <- function(x, y) if (is.null(x)) y else x

resolve_task_refs <- function(x, env) {
    if (is.list(x)) {
        return(lapply(x, resolve_task_refs, env = env))
    }
    if (is.character(x) && length(x) == 1L) {
        if (startsWith(x, "@")) {
            return(eval(parse(text = substring(x, 2L)), envir = env))
        }
        if (startsWith(x, "!!call ")) {
            return(eval(parse(text = substring(x, 7L)), envir = env))
        }
    }
    x
}

.normalize_steps <- function(step = NULL, steps = NULL) {
    x <- steps %||% step
    if (is.null(x)) {
        stop("Task definition is missing both 'step' and 'steps'.")
    }
    if (is.list(x)) {
        x <- unlist(x, use.names = FALSE)
    }
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) {
        stop("Task step(s) must be a non-empty character vector.")
    }
    x
}

.make_task <- function(label, from, steps, final_name, params,
                       store_fast_steps = NULL, params_list = NULL) {
    task <- list(
        label = as.character(label),
        from = as.character(from),
        final_name = as.character(final_name)
    )
    if (length(steps) == 1L) {
        task$step <- as.character(steps[[1]])
    } else {
        task$steps <- as.character(steps)
    }
    if (!is.null(params_list)) {
        task$params_list <- params_list
    } else {
        task$params <- params
    }
    if (!is.null(store_fast_steps)) {
        task$store_fast_steps <- isTRUE(store_fast_steps)
    }
    task
}

.as_named_index <- function(items, section_name) {
    if (!is.list(items) || length(items) == 0L) {
        return(list())
    }
    ids <- vapply(items, function(x) as.character(x$id %||% ""), character(1))
    ids <- trimws(ids)
    if (any(!nzchar(ids))) {
        stop(sprintf("All entries in '%s' must contain a non-empty 'id'.", section_name))
    }
    if (any(duplicated(ids))) {
        dup <- unique(ids[duplicated(ids)])
        stop(sprintf("Duplicate ids in '%s': %s", section_name, paste(dup, collapse = ", ")))
    }
    stats::setNames(items, ids)
}

.is_missing_required_value <- function(x) {
    if (is.null(x)) {
        return(TRUE)
    }
    if (length(x) == 0L) {
        return(TRUE)
    }
    if (length(x) == 1L) {
        if (is.na(x)) {
            return(TRUE)
        }
        if (is.character(x) && !nzchar(trimws(x))) {
            return(TRUE)
        }
    }
    FALSE
}

.check_required_params <- function(params, required_params) {
    required <- required_params %||% character(0)
    required <- as.character(required)
    required <- required[!is.na(required) & nzchar(required)]
    if (length(required) == 0L) {
        return(character(0))
    }
    missing <- required[vapply(required, function(k) {
        !(k %in% names(params)) || .is_missing_required_value(params[[k]])
    }, logical(1))]
    unique(missing)
}

.resolve_params <- function(default_params, method_params, env) {
    merged <- utils::modifyList(default_params %||% list(), method_params %||% list())
    resolve_task_refs(merged, env = env)
}

.read_legacy_tasks <- function(tasks, env) {
    lapply(tasks, function(t) {
        if (is.null(t$from) || (is.null(t$step) && is.null(t$steps)) ||
            (is.null(t$params) && is.null(t$params_list))) {
            stop(
                "Each task must include 'from', one of 'step'/'steps', and one of 'params'/'params_list'."
            )
        }
        if (!is.null(t$params)) {
            t$params <- resolve_task_refs(t$params, env = env)
        }
        if (!is.null(t$params_list)) {
            t$params_list <- resolve_task_refs(t$params_list, env = env)
        }
        if (!is.null(t$steps) && is.list(t$steps)) {
            t$steps <- unlist(t$steps, use.names = FALSE)
        }
        t
    })
}

.append_task <- function(task_list, task_obj) {
    task_list[[length(task_list) + 1L]] <- task_obj
    task_list
}

.prepare_task_env <- function(env) {
    aliased <- new.env(parent = env)

    if (!exists("sample_annotation", envir = env, inherits = TRUE) &&
        exists("sa", envir = env, inherits = TRUE)) {
        assign("sample_annotation", get("sa", envir = env, inherits = TRUE), envir = aliased)
    }

    if (!exists("batch_col", envir = env, inherits = TRUE) &&
        exists("technical_factors", envir = env, inherits = TRUE)) {
        technical_factors <- get("technical_factors", envir = env, inherits = TRUE)
        if (length(technical_factors) > 0L) {
            assign("batch_col", as.character(technical_factors[[1]]), envir = aliased)
        }
    }

    if (!exists("condition_col", envir = env, inherits = TRUE) &&
        exists("biological_factors", envir = env, inherits = TRUE)) {
        biological_factors <- get("biological_factors", envir = env, inherits = TRUE)
        if (length(biological_factors) > 0L) {
            assign("condition_col", as.character(biological_factors[[1]]), envir = aliased)
        }
    }

    if (!exists("assay_prefix", envir = env, inherits = TRUE) &&
        exists("level", envir = env, inherits = TRUE)) {
        assign("assay_prefix", as.character(get("level", envir = env, inherits = TRUE)), envir = aliased)
    }

    aliased
}

.apply_user_params <- function(grid, env) {
    user_params <- grid$user_params %||% list()
    if (!is.list(user_params) || length(user_params) == 0L) {
        return(env)
    }

    keys <- names(user_params)
    if (is.null(keys) || any(!nzchar(keys))) {
        stop("All entries in 'task_grid.user_params' must be named.")
    }
    if (anyDuplicated(keys)) {
        dup <- unique(keys[duplicated(keys)])
        stop(sprintf(
            "Duplicate keys in 'task_grid.user_params': %s",
            paste(dup, collapse = ", ")
        ))
    }

    scoped <- new.env(parent = env)
    for (key in keys) {
        value <- resolve_task_refs(user_params[[key]], env = scoped)
        assign(key, value, envir = scoped)
    }
    scoped
}

.expand_task_grid <- function(grid, env = parent.frame()) {
    if (!is.list(grid)) {
        stop("task_grid must be a YAML mapping.")
    }

    settings <- resolve_task_refs(grid$settings %||% list(), env = env)
    defaults <- grid$defaults %||% list()
    combinations <- resolve_task_refs(grid$combinations %||% list(), env = env)

    assay_prefix <- as.character(settings$assay_prefix %||% "feature")
    input_assay <- as.character(settings$input_assay %||% paste0(assay_prefix, "::log2_on_raw"))
    skip_invalid <- isTRUE(settings$skip_invalid)

    profiles_idx <- .as_named_index(grid$profiles %||% list(), "task_grid.profiles")
    imputers_idx <- .as_named_index(grid$imputation_methods %||% list(), "task_grid.imputation_methods")
    corrections_idx <- .as_named_index(grid$correction_methods %||% list(), "task_grid.correction_methods")

    all_profile_ids <- names(profiles_idx)
    all_imputer_ids <- names(imputers_idx)
    all_correction_ids <- names(corrections_idx)

    profiles_for_direct <- as.character(combinations$profiles_for_direct_correction %||% all_profile_ids)
    profiles_for_imputation <- as.character(combinations$profiles_for_imputation %||% all_profile_ids)
    direct_corrections <- as.character(combinations$direct_corrections %||% all_correction_ids)
    imputation_ids <- as.character(combinations$imputations %||% all_imputer_ids)
    post_imputation_corrections <- as.character(combinations$post_imputation_corrections %||% all_correction_ids)

    validate_ids <- function(ids, known, section_name) {
        missing <- setdiff(unique(ids), known)
        if (length(missing) > 0L) {
            stop(sprintf(
                "Unknown ids in %s: %s",
                section_name,
                paste(missing, collapse = ", ")
            ))
        }
    }
    validate_ids(profiles_for_direct, all_profile_ids, "task_grid.combinations.profiles_for_direct_correction")
    validate_ids(profiles_for_imputation, all_profile_ids, "task_grid.combinations.profiles_for_imputation")
    validate_ids(direct_corrections, all_correction_ids, "task_grid.combinations.direct_corrections")
    validate_ids(imputation_ids, all_imputer_ids, "task_grid.combinations.imputations")
    validate_ids(post_imputation_corrections, all_correction_ids, "task_grid.combinations.post_imputation_corrections")

    tasks <- list()
    profile_assays <- list()

    preprocess_default_params <- defaults$preprocess %||% list()
    imputation_default_params <- defaults$imputation %||% list()
    correction_default_params <- defaults$correction %||% list()
    default_store_fast_steps <- defaults$store_fast_steps

    selected_profile_ids <- unique(c(profiles_for_direct, profiles_for_imputation))
    for (profile_id in selected_profile_ids) {
        profile <- profiles_idx[[profile_id]]
        has_profile_steps <- !is.null(profile$step) || !is.null(profile$steps)
        profile_from <- as.character(profile$from %||% input_assay)
        # Step-less profiles are pass-through anchors that point to an existing assay.
        profile_final <- if (!is.null(profile$final_name)) {
            as.character(profile$final_name)
        } else if (has_profile_steps) {
            as.character(paste0(assay_prefix, "::", profile_id))
        } else {
            profile_from
        }

        if (has_profile_steps) {
            profile_steps <- .normalize_steps(step = profile$step, steps = profile$steps)
            profile_label <- as.character(profile$label %||% profile_id)
            profile_params <- .resolve_params(preprocess_default_params, profile$params, env = env)
            profile_store <- profile$store_fast_steps %||% default_store_fast_steps

            tasks <- .append_task(
                tasks,
                .make_task(
                    label = profile_label,
                    from = profile_from,
                    steps = profile_steps,
                    final_name = profile_final,
                    params = profile_params,
                    store_fast_steps = profile_store
                )
            )
        }

        profile_assays[[profile_id]] <- profile_final
    }

    # Direct correction combinations
    for (profile_id in profiles_for_direct) {
        from_assay <- profile_assays[[profile_id]]
        for (corr_id in direct_corrections) {
            corr <- corrections_idx[[corr_id]]
            corr_steps <- .normalize_steps(step = corr$step, steps = corr$steps)
            corr_params <- .resolve_params(correction_default_params, corr$params, env = env)
            missing_req <- .check_required_params(corr_params, corr$required_params)
            if (length(missing_req) > 0L) {
                msg <- sprintf(
                    "Correction '%s' skipped for profile '%s' (missing required params: %s).",
                    corr_id, profile_id, paste(missing_req, collapse = ", ")
                )
                if (skip_invalid) {
                    message(msg)
                    next
                }
                stop(msg)
            }

            corr_label <- as.character(corr$label %||% paste(corr_id, profile_id, sep = "_"))
            corr_final <- as.character(corr$final_name %||% paste0(
                assay_prefix, "::", corr_id, "_on_", profile_id
            ))
            corr_from <- as.character(corr$from %||% from_assay)

            tasks <- .append_task(
                tasks,
                .make_task(
                    label = corr_label,
                    from = corr_from,
                    steps = corr_steps,
                    final_name = corr_final,
                    params = corr_params
                )
            )
        }
    }

    # Imputation combinations (and post-imputation correction combinations)
    include_imputation_outputs <- isTRUE(combinations$include_imputation_outputs %||% TRUE)
    run_post_imputation <- length(post_imputation_corrections) > 0L
    if (length(imputation_ids) > 0L && (include_imputation_outputs || run_post_imputation)) {
        imputation_assays <- list()
        for (profile_id in profiles_for_imputation) {
            from_assay <- profile_assays[[profile_id]]
            for (imp_id in imputation_ids) {
                imp <- imputers_idx[[imp_id]]
                imp_steps <- .normalize_steps(step = imp$step, steps = imp$steps)
                imp_params <- .resolve_params(imputation_default_params, imp$params, env = env)
                missing_req <- .check_required_params(imp_params, imp$required_params)
                if (length(missing_req) > 0L) {
                    msg <- sprintf(
                        "Imputation '%s' skipped for profile '%s' (missing required params: %s).",
                        imp_id, profile_id, paste(missing_req, collapse = ", ")
                    )
                    if (skip_invalid) {
                        message(msg)
                        next
                    }
                    stop(msg)
                }

                imp_label <- as.character(imp$label %||% paste(imp_id, profile_id, sep = "_"))
                imp_final <- as.character(imp$final_name %||% paste0(
                    assay_prefix, "::", imp_id, "_on_", profile_id
                ))
                imp_from <- as.character(imp$from %||% from_assay)

                tasks <- .append_task(
                    tasks,
                    .make_task(
                        label = imp_label,
                        from = imp_from,
                        steps = imp_steps,
                        final_name = imp_final,
                        params = imp_params
                    )
                )
                key <- paste(profile_id, imp_id, sep = "::")
                imputation_assays[[key]] <- imp_final
            }
        }

        for (key in names(imputation_assays)) {
            parts <- strsplit(key, "::", fixed = TRUE)[[1]]
            profile_id <- parts[[1]]
            imp_id <- parts[[2]]
            imp_assay <- imputation_assays[[key]]

            for (corr_id in post_imputation_corrections) {
                corr <- corrections_idx[[corr_id]]
                corr_steps <- .normalize_steps(step = corr$step, steps = corr$steps)
                corr_params <- .resolve_params(correction_default_params, corr$params, env = env)
                missing_req <- .check_required_params(corr_params, corr$required_params)
                if (length(missing_req) > 0L) {
                    msg <- sprintf(
                        "Post-imputation correction '%s' skipped for '%s' (missing required params: %s).",
                        corr_id, key, paste(missing_req, collapse = ", ")
                    )
                    if (skip_invalid) {
                        message(msg)
                        next
                    }
                    stop(msg)
                }

                corr_label <- as.character(corr$post_label %||% paste(corr_id, imp_id, profile_id, sep = "_"))
                corr_final <- as.character(corr$post_final_name %||% paste0(
                    assay_prefix, "::", corr_id, "_on_", imp_id, "_on_", profile_id
                ))

                tasks <- .append_task(
                    tasks,
                    .make_task(
                        label = corr_label,
                        from = imp_assay,
                        steps = corr_steps,
                        final_name = corr_final,
                        params = corr_params
                    )
                )
            }
        }
    }

    # Deduplicate by final assay name while preserving order.
    dedup <- list()
    seen_final_names <- character()
    for (task in tasks) {
        final_name <- task$final_name %||% ""
        if (nzchar(final_name) && (final_name %in% seen_final_names)) {
            next
        }
        dedup[[length(dedup) + 1L]] <- task
        if (nzchar(final_name)) {
            seen_final_names <- c(seen_final_names, final_name)
        }
    }
    dedup
}

read_pb_tasks_yaml <- function(path, env = parent.frame()) {
    if (!file.exists(path)) {
        stop(sprintf("Task YAML file does not exist: %s", path))
    }
    if (!requireNamespace("yaml", quietly = TRUE)) {
        stop("Package 'yaml' is required to read correction_tasks_yaml.")
    }

    env <- .prepare_task_env(env)
    y <- yaml::read_yaml(path)
    if (!is.null(y$tasks)) {
        return(.read_legacy_tasks(y$tasks, env = env))
    }
    if (!is.null(y$task_grid)) {
        env_with_user_params <- .apply_user_params(y$task_grid, env = env)
        return(.expand_task_grid(y$task_grid, env = env_with_user_params))
    }
    stop("YAML must include either a top-level 'tasks:' list or a 'task_grid:' specification.")
}

filter_pb_tasks <- function(tasks, labels = NULL) {
    if (is.null(labels) || length(labels) == 0L) {
        return(tasks)
    }
    labels <- unique(as.character(labels))
    keep <- vapply(tasks, function(t) {
        !is.null(t$label) && (t$label %in% labels)
    }, logical(1))
    tasks[keep]
}

build_simple_correction_tasks <- function(
  correction_methods,
  assay_prefix,
  batch_col,
  correction_covariates = character(0),
  from_assay = paste0(assay_prefix, "::log2_on_raw")
) {
    methods <- unique(as.character(correction_methods))
    methods <- methods[!is.na(methods) & nzchar(methods)]
    if (length(methods) == 0L) {
        stop("correction_methods must contain at least one non-empty method name.")
    }

    lapply(methods, function(method) {
        params <- list(batch_col = batch_col)
        if (length(correction_covariates) > 0L) {
            params$covariates_cols <- correction_covariates
        }
        list(
            label = method,
            from = from_assay,
            step = method,
            final_name = paste0(assay_prefix, "::", method, "_on_log2_on_raw"),
            params = params
        )
    })
}

sanitize_label_for_path <- function(x) {
    if (is.null(x) || !nzchar(x)) {
        return("unnamed_method")
    }
    gsub("[^A-Za-z0-9._-]+", "_", x)
}

run_pb_tasks <- function(pbf, tasks, log_fn = message) {
    if (length(tasks) == 0L) {
        return(list(
            pbf = pbf,
            task_results = data.frame(),
            corrected_assays = character(0)
        ))
    }

    task_rows <- vector("list", length(tasks))

    for (i in seq_along(tasks)) {
        t <- tasks[[i]]
        steps <- t$steps %||% t$step
        if (is.list(steps)) {
            steps <- unlist(steps, use.names = FALSE)
        }
        if (!is.character(steps) || length(steps) < 1L) {
            stop("Task step(s) must be a non-empty character vector.")
        }

        params_list <- t$params_list
        if (is.null(params_list)) {
            p <- t$params
            if (length(steps) == 1L) {
                params_list <- list(p)
            } else if (is.list(p) && length(p) == length(steps) &&
                all(vapply(p, is.list, logical(1)))) {
                params_list <- p
            } else if (!is.null(p)) {
                params_list <- rep(list(p), length(steps))
            } else {
                stop("Task requires 'params_list' or reusable 'params'.")
            }
        }

        label <- t$label %||% paste(steps, collapse = "_")
        step_text <- paste(steps, collapse = " -> ")
        log_fn(sprintf("(%d/%d) Running task '%s' [%s]", i, length(tasks), label, step_text))

        transform_args <- list(
            pbf,
            from = t$from,
            steps = steps,
            store_fast_steps = isTRUE(t$store_fast_steps),
            params_list = params_list
        )
        if (!is.null(t$final_name)) {
            transform_args$final_name <- t$final_name
        }

        task_error <- NULL
        updated <- tryCatch(
            do.call(proBatch::pb_transform, transform_args),
            error = function(e) {
                task_error <<- conditionMessage(e)
                pbf
            }
        )

        if (is.null(task_error)) {
            pbf <- updated
        }

        final_name <- t$final_name
        if (is.null(final_name) && is.null(task_error)) {
            final_name <- proBatch::pb_current_assay(pbf)
        }
        assay_created <- !is.null(final_name) && (final_name %in% names(pbf))
        success <- is.null(task_error) && assay_created

        if (!is.null(task_error)) {
            log_fn(sprintf("Task '%s' failed: %s", label, task_error))
        } else if (!assay_created) {
            log_fn(sprintf("Task '%s' finished but assay '%s' was not found.", label, final_name))
        }

        task_rows[[i]] <- data.frame(
            task_index = i,
            label = as.character(label),
            from = as.character(t$from),
            steps = as.character(step_text),
            final_name = as.character(final_name %||% NA_character_),
            success = success,
            message = if (is.null(task_error)) {
                if (assay_created) "ok" else "assay_not_found"
            } else {
                task_error
            },
            stringsAsFactors = FALSE
        )
    }

    task_results <- do.call(rbind, task_rows)
    corrected_assays <- unique(task_results$final_name[task_results$success])
    corrected_assays <- corrected_assays[!is.na(corrected_assays) & nzchar(corrected_assays)]

    list(
        pbf = pbf,
        task_results = task_results,
        corrected_assays = corrected_assays
    )
}

# ==============================================================================
# Workflow-level reusable helpers
# ==============================================================================

`%+%` <- function(x, y) paste0(x, y)

log_msg <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

apply_workflow_params <- function(workflow_params, env = parent.frame(), log_fn = log_msg) {
    if (is.null(workflow_params)) {
        return(invisible(NULL))
    }
    if (!is.list(workflow_params)) {
        stop("'workflow_params' must be a named list.")
    }
    nms <- names(workflow_params)
    if (is.null(nms) || any(!nzchar(nms))) {
        stop("'workflow_params' must be a named list with non-empty keys.")
    }

    unknown <- setdiff(nms, ls(envir = env, all.names = TRUE))
    if (length(unknown) > 0L) {
        log_fn(sprintf(
            "workflow_params contains unknown key(s): %s",
            paste(unknown, collapse = ", ")
        ), level = "WARN")
    }

    for (key in nms) {
        assign(key, workflow_params[[key]], envir = env)
    }

    log_fn(sprintf("Applied %d workflow parameter override(s).", length(nms)))
    invisible(NULL)
}

safe_try <- function(expr, error_msg = "Error occurred", return_on_error = NULL) {
    result <- try(expr, silent = TRUE)
    if (inherits(result, "try-error")) {
        log_msg(sprintf("%s: %s", error_msg, as.character(result)), level = "WARN")
        return(return_on_error)
    }
    result
}

save_plot_helper <- function(
  plot_obj,
  filename,
  subdir = "plots",
  width = plot_width,
  height = plot_height,
  format = plot_format,
  dpi = plot_dpi,
  n_plots = NULL,
  n_cols = NULL
) {
    if (!save_plots) {
        return(invisible(NULL))
    }

    filepath <- file.path(subdir, sprintf("%s.%s", filename, format))

    # If multiple plots are arranged, scale device size to match the facet grid.
    if (!is.null(n_plots)) {
        n_plots_num <- if (is.numeric(n_plots)) n_plots else length(n_plots)
        if (n_plots_num > 1) {
            # Keep default grid logic unless caller explicitly provides n_cols.
            n_cols_use <- n_cols
            if (is.null(n_cols_use)) {
                n_cols_use <- if (n_plots_num > 8) 8 else if (n_plots_num > 5) 5 else 2
            }
            n_rows <- ceiling(n_plots_num / n_cols_use)
            width <- round(width * n_cols_use / 2, 0)
            height <- round(height * n_rows, 0)
        }
    }

    safe_try(
        {
            if (format == "pdf") {
                pdf(filepath, width = width, height = height)
            } else if (format == "png") {
                png(filepath, width = width, height = height, units = "in", res = dpi)
            } else if (format == "svg") {
                svg(filepath, width = width, height = height)
            }
            print(plot_obj)
            dev.off()
            log_msg(sprintf("Saved plot: %s", filepath))
        },
        error_msg = sprintf("Failed to save plot: %s", filepath)
    )
}

save_metric_helper <- function(data_obj, filename, subdir = "metrics") {
    if (!save_metrics) {
        return(invisible(NULL))
    }

    filepath <- file.path(subdir, sprintf("%s.csv", filename))

    safe_try(
        {
            if (is.data.frame(data_obj) || is.matrix(data_obj)) {
                write.csv(data_obj, filepath, row.names = TRUE)
                log_msg(sprintf("Saved metric: %s", filepath))
            } else {
                safe_try(
                    {
                        metric_flatten <- pb_flatten_design_check(data_obj)
                        write.csv(metric_flatten, filepath, row.names = FALSE)
                        log_msg(sprintf("Saved metric (flattened): %s", filepath))
                    },
                    error_msg = sprintf("Failed to save flattened metric: %s", filepath)
                )
            }
        },
        error_msg = sprintf("Failed to save metric: %s", filepath)
    )
}

safe_classification_metrics <- function(
  x,
  sample_annotation,
  known_col = batch_col,
  fill_the_missing = fill_missing
) {
    if (is.null(sample_annotation) || !all(known_col %in% names(sample_annotation))) {
        return(NULL)
    }

    for (col in known_col) {
        n_levels <- length(unique(sample_annotation[[col]][!is.na(sample_annotation[[col]])]))
        if (n_levels < 2L) {
            log_msg(sprintf(
                "Not enough levels in column '%s' for classification metrics. Skipping.",
                col
            ), level = "WARN")
            known_col <- setdiff(known_col, col)
        }
    }
    if (length(known_col) == 0) {
        log_msg("No valid known columns for classification metrics. Skipping.", level = "WARN")
        return(NULL)
    }

    safe_try(
        proBatch::calculate_classification_metrics(
            x,
            known_col = known_col,
            fill_the_missing = fill_the_missing
        ),
        error_msg = "Classification metrics calculation failed",
        return_on_error = NULL
    )
}

pb_flatten_design_check <- function(x, sep = ", ") {
    `%||%` <- function(a, b) if (is.null(a)) b else a
    empty <- data.frame(section = character(), key = character(), value = character(), stringsAsFactors = FALSE)
    if (is.null(x)) {
        return(empty)
    }

    if (!is.list(x)) {
        return(data.frame(
            section = "value",
            key = "",
            value = paste(x, collapse = sep),
            stringsAsFactors = FALSE
        ))
    }

    summ <- x$summary %||% x[setdiff(names(x), c("errors", "warnings", "summary"))]
    summ <- if (length(summ)) {
        rapply(summ, function(v) if (length(v)) paste(v, collapse = sep) else NA_character_, how = "replace")
    } else {
        list()
    }
    flat <- if (length(summ)) unlist(summ, recursive = TRUE, use.names = TRUE) else character(0)

    nm <- names(flat)
    if (is.null(nm)) nm <- rep("", length(flat))
    parts <- strsplit(nm, "\\.", perl = TRUE)
    sec <- vapply(parts, function(p) if (length(p) && nzchar(p[1])) p[1] else "summary", character(1))
    key <- vapply(parts, function(p) if (length(p) > 1) paste(p[-1], collapse = ".") else "", character(1))

    errs <- as.character(unlist(x$errors %||% character(0), use.names = FALSE))
    warns <- as.character(unlist(x$warnings %||% character(0), use.names = FALSE))

    out <- do.call(rbind, Filter(Negate(is.null), list(
        if (length(errs)) {
            data.frame(section = "errors", key = as.character(seq_along(errs)), value = errs, stringsAsFactors = FALSE)
        } else {
            NULL
        },
        if (length(warns)) {
            data.frame(section = "warnings", key = as.character(seq_along(warns)), value = warns, stringsAsFactors = FALSE)
        } else {
            NULL
        },
        if (length(flat)) {
            data.frame(section = sec, key = key, value = as.character(unname(flat)), stringsAsFactors = FALSE)
        } else {
            NULL
        }
    )))
    if (is.null(out)) empty else out
}

first_num <- function(x) {
    if (is.null(x) || length(x) == 0) {
        return(NA_real_)
    }
    x <- suppressWarnings(as.numeric(x))
    if (!length(x) || all(is.na(x))) {
        return(NA_real_)
    }
    x[which(!is.na(x))[1]]
}

pvca_share <- function(pvca_df, patterns = character(), category = NULL) {
    if (is.null(pvca_df) || !is.data.frame(pvca_df)) {
        return(NA_real_)
    }

    value_col <- if ("weights" %in% names(pvca_df)) {
        "weights"
    } else if ("Weighted.average.proportion.variance" %in% names(pvca_df)) {
        "Weighted.average.proportion.variance"
    } else {
        return(NA_real_)
    }

    matched <- pvca_df
    if (!is.null(category) && "category" %in% names(pvca_df)) {
        matched <- pvca_df[tolower(pvca_df$category) == tolower(category), , drop = FALSE]
    } else if (!is.null(patterns) && length(patterns) > 0 && "label" %in% names(pvca_df)) {
        pattern_vec <- patterns[!is.na(patterns) & nzchar(patterns)]
        if (length(pattern_vec) == 0) {
            return(NA_real_)
        }
        matched <- pvca_df[grepl(
            paste(pattern_vec, collapse = "|"),
            pvca_df$label,
            ignore.case = TRUE
        ), , drop = FALSE]
    } else {
        return(NA_real_)
    }

    if (nrow(matched) == 0) {
        return(NA_real_)
    }
    sum(matched[[value_col]], na.rm = TRUE)
}

median_corr_by_group <- function(corr_df, group_label) {
    if (is.null(corr_df) || !is.data.frame(corr_df) || !"correlation" %in% names(corr_df)) {
        return(NA_real_)
    }

    idx <- rep(FALSE, nrow(corr_df))
    if (identical(group_label, "within_replicate")) {
        if ("replicate" %in% names(corr_df)) {
            replicate_flag <- suppressWarnings(as.logical(corr_df$replicate))
            idx <- !is.na(replicate_flag) & replicate_flag
        } else if ("batch_replicate" %in% names(corr_df)) {
            idx <- grepl("same_biospecimen|within_replicate", corr_df$batch_replicate, ignore.case = TRUE)
        }
    } else if (identical(group_label, "within_batch")) {
        if ("batch_the_same" %in% names(corr_df)) {
            batch_flag <- suppressWarnings(as.logical(corr_df$batch_the_same))
            idx <- !is.na(batch_flag) & batch_flag
        } else if ("batch_replicate" %in% names(corr_df)) {
            idx <- grepl("same_batch|within_batch", corr_df$batch_replicate, ignore.case = TRUE)
        }
    }

    if (!any(idx, na.rm = TRUE)) {
        return(NA_real_)
    }
    median(corr_df$correlation[idx], na.rm = TRUE)
}

count_outliers <- function(x) {
    if (is.null(x)) {
        return(NA_integer_)
    }
    if (is.data.frame(x)) {
        outlier_col <- c("is_outlier", "outlier")
        outlier_col <- outlier_col[outlier_col %in% names(x)]
        if (length(outlier_col) > 0) {
            outlier_flag <- suppressWarnings(as.logical(x[[outlier_col[1]]]))
            return(sum(outlier_flag, na.rm = TRUE))
        }
    }
    if (is.list(x) && "outliers" %in% names(x)) {
        return(length(x$outliers))
    }
    if (is.logical(x)) {
        return(sum(x, na.rm = TRUE))
    }
    NA_integer_
}

extract_class_metric <- function(class_df, metric_cols, known_col_target = batch_col) {
    if (is.null(class_df) || !is.data.frame(class_df)) {
        return(NA_real_)
    }
    class_subset <- class_df
    if (!is.null(known_col_target) && "known_col" %in% names(class_subset) &&
        known_col_target %in% class_subset$known_col) {
        class_subset <- class_subset[class_subset$known_col == known_col_target, , drop = FALSE]
    }

    metric_cols <- metric_cols[metric_cols %in% names(class_subset)]
    if (length(metric_cols) == 0) {
        return(NA_real_)
    }
    first_num(class_subset[[metric_cols[1]]])
}

median_feature_cv <- function(cv_df) {
    if (is.null(cv_df) || !is.data.frame(cv_df)) {
        return(NA_real_)
    }
    cv_col <- c("CV_replicate", "CV_total", "CV_perBatch")
    cv_col <- cv_col[cv_col %in% names(cv_df)]
    if (length(cv_col) == 0) {
        return(NA_real_)
    }
    cv_values <- suppressWarnings(as.numeric(cv_df[[cv_col[1]]]))
    if (!length(cv_values) || all(is.na(cv_values))) {
        return(NA_real_)
    }
    median(cv_values, na.rm = TRUE)
}
