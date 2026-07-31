# Provider-neutral transformation registry ---------------------------------

.pb_step_records <- new.env(parent = emptyenv())
.pb_step_aliases <- new.env(parent = emptyenv())

.pb_registry_character_scalar <- function(value, argument) {
    if (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(value)) {
        stop(
            "`", argument,
            "` must be one non-missing, non-empty character value.",
            call. = FALSE
        )
    }
    value
}

.pb_registry_character_vector <- function(value, argument) {
    if (!is.character(value) || anyNA(value) || any(!nzchar(value))) {
        stop(
            "`", argument,
            "` must be a character vector of non-missing, non-empty values.",
            call. = FALSE
        )
    }
    if (anyDuplicated(value)) {
        stop("`", argument, "` must not contain duplicates.", call. = FALSE)
    }
    value
}

.pb_step_environment_name <- function(fun) {
    fun_environment <- environment(fun)
    if (is.null(fun_environment)) {
        return("base")
    }
    environmentName(fun_environment)
}

.pb_infer_step_package <- function(fun) {
    fun_environment <- environment(fun)
    if (is.null(fun_environment)) {
        return("base")
    }

    environment_name <- environmentName(fun_environment)
    if (!nzchar(environment_name) ||
        environment_name %in% c("R_GlobalEnv", "Autoloads")) {
        return(NA_character_)
    }

    package_name <- sub("^(namespace|package):", "", environment_name)
    if (!nzchar(package_name)) {
        return(NA_character_)
    }
    package_name
}

.pb_normalize_step_package <- function(package, fun) {
    if (is.null(package)) {
        return(.pb_infer_step_package(fun))
    }
    .pb_registry_character_scalar(package, "package")
}

.pb_step_record_identical <- function(left, right) {
    fields <- c(
        "name", "fun", "package", "kind", "requires", "aliases", "label"
    )
    all(vapply(
        fields,
        function(field) identical(left[[field]], right[[field]]),
        logical(1)
    ))
}

.pb_lookup_step_record <- function(name) {
    if (exists(name, envir = .pb_step_records, inherits = FALSE)) {
        return(get(name, envir = .pb_step_records, inherits = FALSE))
    }
    if (!exists(name, envir = .pb_step_aliases, inherits = FALSE)) {
        return(NULL)
    }
    canonical_name <- get(name, envir = .pb_step_aliases, inherits = FALSE)
    get(canonical_name, envir = .pb_step_records, inherits = FALSE)
}

.pb_requirement_available <- function(package) {
    length(find.package(package, quiet = TRUE)) == 1L
}

.pb_step_availability <- function(record) {
    provider_missing <- !is.na(record$package) &&
        !record$package %in% loadedNamespaces()
    missing_requirements <- record$requires[
        !vapply(record$requires, .pb_requirement_available, logical(1))
    ]

    list(
        available = !provider_missing && !length(missing_requirements),
        provider_missing = provider_missing,
        missing_requirements = missing_requirements
    )
}

.pb_step_provider_label <- function(package) {
    if (is.na(package)) {
        return("<ordinary registration>")
    }
    package
}

.pb_expected_step_package <- function(package) {
    if (is.null(package)) {
        return(NULL)
    }
    if (!is.character(package) || length(package) != 1L ||
        (!is.na(package) && !nzchar(package))) {
        stop(
            "`package` must be NULL, NA, or one non-empty character value.",
            call. = FALSE
        )
    }
    package
}

.pb_stop_step_unavailable <- function(record, availability) {
    reasons <- character()
    if (isTRUE(availability$provider_missing)) {
        reasons <- c(
            reasons,
            paste0(
                "provider namespace '", record$package, "' is not loaded"
            )
        )
    }
    if (length(availability$missing_requirements)) {
        reasons <- c(
            reasons,
            paste0(
                "requires unavailable package",
                if (length(availability$missing_requirements) ==
                    1L) {
                    ""
                } else {
                    "s"
                },
                ": ",
                paste(
                    shQuote(availability$missing_requirements),
                    collapse = ", "
                )
            )
        )
    }

    stop(
        "Registered step '", record$name, "' is unavailable: ",
        paste(reasons, collapse = "; "), ". ",
        "Load the provider and make its requirements available, then register ",
        "the provider again if needed.",
        call. = FALSE
    )
}

#' Register a provider-neutral transformation step
#'
#' Registers one canonical step and optional aliases for use by
#' [pb_transform()], [pb_eval()], and [pb_apply_matrix_method()]. Registration
#' records metadata only: it does not install, load, or probe a provider or any
#' optional engine. Repeating an identical registration is a no-op. A different
#' registration that collides with an existing canonical name or alias is
#' rejected; refresh a provider with [pb_unregister_steps()]
#' followed by registration.
#'
#' @param name One non-missing, non-empty canonical step name.
#' @param fun Function implementing the step. It receives a numeric matrix as
#'   its first argument and may declare a `sample_annotation` argument. It must
#'   return a numeric matrix or a [pb_step_result()].
#' @param package Optional provider package name. When omitted, the function
#'   environment is used to infer a provider without loading a package.
#'   Functions without a package namespace are ordinary user registrations.
#' @param kind One non-missing, non-empty step-kind label.
#' @param requires Character vector naming packages required when the step is
#'   invoked.
#' @param aliases Character vector of alternative names for `name`.
#' @param label One non-missing, non-empty display label.
#'
#' @return `TRUE`, invisibly.
#' @example inst/examples/ProBatchFeatures-basic.R
#' @export
#' @md
pb_register_step <- function(
  name,
  fun,
  package = NULL,
  kind = "transform",
  requires = character(),
  aliases = character(),
  label = name
) {
    name <- .pb_registry_character_scalar(name, "name")
    if (!is.function(fun)) {
        stop("`fun` must be a function.", call. = FALSE)
    }
    package <- .pb_normalize_step_package(package, fun)
    kind <- .pb_registry_character_scalar(kind, "kind")
    requires <- .pb_registry_character_vector(requires, "requires")
    aliases <- .pb_registry_character_vector(aliases, "aliases")
    label <- .pb_registry_character_scalar(label, "label")

    if (name %in% aliases) {
        stop("`aliases` must not contain the canonical `name`.", call. = FALSE)
    }

    record <- list(
        name = name,
        fun = fun,
        package = package,
        kind = kind,
        requires = requires,
        aliases = aliases,
        label = label
    )

    if (exists(name, envir = .pb_step_records, inherits = FALSE)) {
        existing <- get(name, envir = .pb_step_records, inherits = FALSE)
        if (.pb_step_record_identical(existing, record)) {
            return(invisible(TRUE))
        }
        stop(
            "Canonical step '", name,
            "' is already registered with different metadata or function. ",
            "Unregister its provider before registering a replacement.",
            call. = FALSE
        )
    }

    if (exists(name, envir = .pb_step_aliases, inherits = FALSE)) {
        owner <- get(name, envir = .pb_step_aliases, inherits = FALSE)
        stop(
            "Canonical step '", name, "' collides with alias '", name,
            "' owned by canonical step '", owner, "'.",
            call. = FALSE
        )
    }

    for (alias in aliases) {
        if (exists(alias, envir = .pb_step_records, inherits = FALSE)) {
            stop(
                "Alias '", alias, "' collides with canonical step '", alias,
                "'.",
                call. = FALSE
            )
        }
        if (exists(alias, envir = .pb_step_aliases, inherits = FALSE)) {
            owner <- get(alias, envir = .pb_step_aliases, inherits = FALSE)
            stop(
                "Alias '", alias, "' is already owned by canonical step '",
                owner, "'.",
                call. = FALSE
            )
        }
    }

    assign(name, record, envir = .pb_step_records)
    for (alias in aliases) {
        assign(alias, name, envir = .pb_step_aliases)
    }
    invisible(TRUE)
}

#' Unregister all transformation steps from one provider
#'
#' Removes every canonical registration and alias owned by `package`. Calling
#' the function when the provider owns no registrations is a no-op.
#'
#' @param package One non-missing, non-empty provider package name.
#'
#' @return The removed canonical names, invisibly.
#' @examples
#' local({
#'     provider <- "proBatchExampleProvider"
#'     on.exit(pb_unregister_steps(provider), add = TRUE)
#'     pb_register_step(
#'         "example_unregister_identity",
#'         function(x) x,
#'         package = provider,
#'         aliases = "example_unregister_no_change"
#'     )
#'     alias_registered <- pb_has_step("example_unregister_no_change")
#'     removed <- pb_unregister_steps(provider)
#'     c(
#'         alias_registered = alias_registered,
#'         alias_removed = !pb_has_step("example_unregister_no_change"),
#'         canonical_removed = identical(
#'             removed,
#'             "example_unregister_identity"
#'         )
#'     )
#' })
#' @export
#' @md
pb_unregister_steps <- function(package) {
    package <- .pb_registry_character_scalar(package, "package")
    canonical_names <- ls(envir = .pb_step_records, all.names = TRUE)
    remove_names <- canonical_names[vapply(
        canonical_names,
        function(name) {
            record <- get(name, envir = .pb_step_records, inherits = FALSE)
            !is.na(record$package) && identical(record$package, package)
        },
        logical(1)
    )]
    remove_names <- sort(remove_names)

    if (!length(remove_names)) {
        return(invisible(character()))
    }

    aliases <- unlist(
        lapply(
            remove_names,
            function(name) {
                get(name, envir = .pb_step_records, inherits = FALSE)$aliases
            }
        ),
        use.names = FALSE
    )
    if (length(aliases)) {
        rm(list = aliases, envir = .pb_step_aliases)
    }
    rm(list = remove_names, envir = .pb_step_records)
    invisible(remove_names)
}

#' List registered transformation steps
#'
#' Returns canonical step names. Pattern matching considers canonical names,
#' aliases, and display labels. Availability filtering reports current state
#' without changing registrations.
#'
#' @param pattern Optional regular expression used to filter canonical names,
#'   aliases, and labels.
#' @param details Logical scalar. If `FALSE`, return canonical names. If `TRUE`,
#'   return an `S4Vectors::DataFrame` with provider metadata and
#'   legacy inspection columns.
#' @param available Optional logical scalar. Use `TRUE` for currently available
#'   registrations, `FALSE` for unavailable registrations, or `NULL` for both.
#'
#' @return A character vector of canonical names, or an
#'   `S4Vectors::DataFrame` when `details = TRUE`. Detailed output contains
#'   `step`, `pkg`, `env`, and `n_formals` for compatibility, plus `name`,
#'   `package`, `kind`, `label`, `requires`, `aliases`, and `available`.
#' @examples
#' pb_list_steps()
#' pb_list_steps(details = TRUE)
#' pb_list_steps("^l")
#' @export
#' @md
pb_list_steps <- function(pattern = NULL, details = FALSE, available = NULL) {
    if (!is.logical(details) || length(details) != 1L || is.na(details)) {
        stop("`details` must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(pattern)) {
        pattern <- .pb_registry_character_scalar(pattern, "pattern")
    }
    if (!is.null(available) &&
        (!is.logical(available) || length(available) != 1L ||
            is.na(available))) {
        stop("`available` must be NULL, TRUE, or FALSE.", call. = FALSE)
    }

    canonical_names <- sort(ls(envir = .pb_step_records, all.names = TRUE))
    records <- lapply(
        canonical_names,
        get,
        envir = .pb_step_records,
        inherits = FALSE
    )

    if (!is.null(pattern)) {
        keep <- vapply(
            records,
            function(record) {
                any(grepl(
                    pattern,
                    c(record$name, record$aliases, record$label)
                ))
            },
            logical(1)
        )
        canonical_names <- canonical_names[keep]
        records <- records[keep]
    }

    if (!details && is.null(available)) {
        return(canonical_names)
    }

    availability <- vapply(
        records,
        function(record) .pb_step_availability(record)$available,
        logical(1)
    )
    if (!is.null(available)) {
        keep <- availability == available
        canonical_names <- canonical_names[keep]
        records <- records[keep]
        availability <- availability[keep]
    }

    if (!details) {
        return(canonical_names)
    }

    functions <- lapply(records, `[[`, "fun")
    packages <- vapply(records, `[[`, character(1), "package")
    requires <- lapply(records, `[[`, "requires")
    aliases <- lapply(records, `[[`, "aliases")

    S4Vectors::DataFrame(
        step = canonical_names,
        pkg = packages,
        env = vapply(
            functions,
            .pb_step_environment_name,
            character(1)
        ),
        n_formals = vapply(
            functions,
            function(fun) length(formals(fun)),
            integer(1)
        ),
        name = canonical_names,
        package = packages,
        kind = vapply(records, `[[`, character(1), "kind"),
        label = vapply(records, `[[`, character(1), "label"),
        requires = I(requires),
        aliases = I(aliases),
        available = availability,
        row.names = NULL
    )
}

#' Test whether a transformation step is registered or available
#'
#' Canonical names and aliases resolve to the same registration.
#'
#' @param name One canonical step name or alias.
#' @param available Logical scalar. If `FALSE`, test registration only. If
#'   `TRUE`, also require any recorded provider to be loaded and all declared
#'   requirements to be available.
#'
#' @return A logical scalar.
#' @examples
#' pb_has_step("medianNorm")
#' pb_has_step("combat", available = TRUE)
#' @export
#' @md
pb_has_step <- function(name, available = FALSE) {
    if (!is.logical(available) || length(available) != 1L ||
        is.na(available)) {
        stop("`available` must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.character(name) || length(name) != 1L ||
        is.na(name) || !nzchar(name)) {
        return(FALSE)
    }

    record <- .pb_lookup_step_record(name)
    if (is.null(record)) {
        return(FALSE)
    }
    if (!available) {
        return(TRUE)
    }
    .pb_step_availability(record)$available
}

.pb_resolve_step <- function(
  name,
  package = NULL,
  require_available = TRUE
) {
    package <- .pb_expected_step_package(package)
    if (!is.logical(require_available) || length(require_available) != 1L ||
        is.na(require_available)) {
        stop("`require_available` must be TRUE or FALSE.", call. = FALSE)
    }

    if (is.function(name)) {
        inferred_package <- .pb_infer_step_package(name)
        if (!is.null(package) && !identical(inferred_package, package)) {
            stop(
                "Direct function comes from provider '",
                .pb_step_provider_label(inferred_package),
                "', not requested provider '", package, "'.",
                call. = FALSE
            )
        }
        return(list(
            input_name = name,
            name = NA_character_,
            fun = name,
            package = inferred_package,
            kind = "transform",
            requires = character(),
            aliases = character(),
            label = NA_character_,
            available = TRUE,
            registered = FALSE
        ))
    }

    name <- .pb_registry_character_scalar(name, "name")
    record <- .pb_lookup_step_record(name)
    if (is.null(record)) {
        provider_guidance <- if (is.null(package)) {
            "Register it with `pb_register_step()` before invoking it."
        } else if (is.na(package)) {
            paste0(
                "Step '", name,
                "' recorded as an ordinary registration is not registered. ",
                "Register the same function with `pb_register_step()` before ",
                "replay."
            )
        } else {
            paste0(
                "Step '", name, "' recorded from provider '", package,
                "' is not registered. Load that provider and register it ",
                "with `pb_register_step()` before replay."
            )
        }
        if (!is.null(package)) {
            stop(provider_guidance, call. = FALSE)
        }
        stop(
            "Step '", name, "' is not registered. ", provider_guidance,
            call. = FALSE
        )
    }

    if (!is.null(package) && !identical(record$package, package)) {
        expected_provider <- if (is.na(package)) {
            "an ordinary registration"
        } else {
            paste0("recorded provider '", package, "'")
        }
        stop(
            "Step '", name, "' resolves to canonical step '", record$name,
            "' from provider '", .pb_step_provider_label(record$package),
            "', not ", expected_provider, ". ",
            "Load and register the recorded provider before replay.",
            call. = FALSE
        )
    }

    availability <- .pb_step_availability(record)
    if (require_available && !availability$available) {
        .pb_stop_step_unavailable(record, availability)
    }

    c(
        list(input_name = name),
        record,
        list(
            available = availability$available,
            registered = TRUE
        )
    )
}

.pb_get_step_fun <- function(fun_or_name, use_registry = TRUE) {
    if (is.function(fun_or_name)) {
        return(fun_or_name)
    }
    if (isTRUE(use_registry) &&
        is.character(fun_or_name) && length(fun_or_name) == 1L &&
        !is.na(fun_or_name) && nzchar(fun_or_name) &&
        pb_has_step(fun_or_name)) {
        return(.pb_resolve_step(fun_or_name)$fun)
    }

    # Preserve direct exact-name Core and base-function lookup for existing
    # pipelines whose functions predate provider registration.
    core_namespace <- topenv(environment())
    if (is.character(fun_or_name) && length(fun_or_name) == 1L &&
        !is.na(fun_or_name) && nzchar(fun_or_name) &&
        exists(fun_or_name, envir = core_namespace, inherits = FALSE)) {
        return(get(fun_or_name, envir = core_namespace, inherits = FALSE))
    }
    match.fun(fun_or_name)
}

.pb_builtin_log2_step <- function(
  m,
  pseudo = NULL,
  log_base = 2,
  offset = NULL
) {
    log_transform_dm.default(
        m,
        log_base = log_base,
        offset = offset %||% pseudo %||% 1
    )
}

.pb_builtin_log_step <- function(
  m,
  base = exp(1),
  pseudo = NULL,
  log_base = NULL,
  offset = NULL
) {
    log_transform_dm.default(
        m,
        log_base = log_base %||% base,
        offset = offset %||% pseudo %||% 1
    )
}

.pb_builtin_median_norm_step <- function(
  m,
  sample_annotation = NULL,
  sample_id_col = "FullRunName",
  group_col = NULL,
  inside_batch = FALSE,
  fill_the_missing = "keep",
  fill_value = NULL
) {
    missing_policy <- .pb_normalize_missing_policy(
        fill_the_missing,
        fill_value = fill_value,
        argument = "fill_the_missing"
    )
    if (anyNA(m)) {
        m <- handle_missing_values(
            data_matrix = m,
            warning_message = paste(
                "Median normalization: applying requested missing value",
                "handling before centering."
            ),
            fill_the_missing = missing_policy$policy,
            fill_value = missing_policy$fill_value
        )
        if (!nrow(m) || !ncol(m)) {
            stop(
                "No data remaining after handling missing values for ",
                "median normalization",
                call. = FALSE
            )
        }
    }
    normalize_sample_medians_dm(
        data_matrix = m,
        sample_annotation = sample_annotation,
        sample_id_col = sample_id_col,
        group_col = group_col,
        inside_batch = inside_batch
    )
}
