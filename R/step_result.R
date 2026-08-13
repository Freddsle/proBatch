#' Return transformed data with structured method artifacts
#'
#' Provider-neutral transformation methods can use this result type to return
#' transformed data together with artifacts such as model fits, latent factors,
#' or diagnostics. Matrix-oriented Core APIs validate the \code{data} component
#' independently of the artifacts.
#'
#' @param data Transformed data produced by a method.
#' @param artifacts A list of structured artifacts associated with \code{data}.
#'
#' @return An object of class \code{pb_step_result} with \code{data} and
#'   \code{artifacts} components.
#' @export
#' @examples
#' result <- pb_step_result(
#'     matrix(1:4, nrow = 2),
#'     artifacts = list(converged = TRUE)
#' )
#' result$data
#' result$artifacts
pb_step_result <- function(data, artifacts = list()) {
    if (!is.list(artifacts)) {
        stop("`artifacts` must be a list.", call. = FALSE)
    }

    structure(
        list(data = data, artifacts = artifacts),
        class = "pb_step_result"
    )
}

.pb_step_artifacts_key <- "pb_step_artifacts"

.pb_step_result_parts <- function(value) {
    if (!inherits(value, "pb_step_result")) {
        return(list(
            data = value,
            artifacts = list(),
            structured = FALSE
        ))
    }

    required <- c("data", "artifacts")
    valid_names <- is.list(value) && identical(names(value), required)
    if (!valid_names || !is.list(value$artifacts)) {
        stop(
            "Invalid `pb_step_result`: expected one `data` component and ",
            "one list-valued `artifacts` component.",
            call. = FALSE
        )
    }

    list(
        data = value$data,
        artifacts = value$artifacts,
        structured = TRUE
    )
}

.pb_step_result_matrix_parts <- function(value, context) {
    result <- .pb_step_result_parts(value)
    valid_matrix <- is.matrix(result$data) && is.numeric(result$data)
    if (!valid_matrix) {
        if (result$structured) {
            stop(
                context,
                " must contain a numeric matrix in `data`.",
                call. = FALSE
            )
        }
        stop(context, " must be a numeric matrix.", call. = FALSE)
    }
    result
}

.pb_attach_step_artifacts <- function(assay, artifacts) {
    if (!is.list(artifacts)) {
        stop("`artifacts` must be a list.", call. = FALSE)
    }
    assay_metadata <- S4Vectors::metadata(assay)
    assay_metadata[[.pb_step_artifacts_key]] <- artifacts
    S4Vectors::metadata(assay) <- assay_metadata
    assay
}

.pb_read_step_artifact_metadata <- function(assay, context) {
    assay_metadata <- S4Vectors::metadata(assay)
    present <- .pb_step_artifacts_key %in% names(assay_metadata)
    if (!present) {
        return(list(present = FALSE, artifacts = list()))
    }

    artifacts <- assay_metadata[[.pb_step_artifacts_key]]
    if (!is.list(artifacts)) {
        stop(
            context,
            " has malformed structured artifact metadata; expected a list.",
            call. = FALSE
        )
    }
    list(present = TRUE, artifacts = artifacts)
}

.pb_assert_step_artifacts_match <- function(assay, expected, target) {
    stored <- .pb_read_step_artifact_metadata(
        assay,
        paste0("Result target '", target, "'")
    )
    if (!stored$present || !identical(stored$artifacts, expected)) {
        stop(
            "Result target '",
            target,
            "' already exists with conflicting structured artifacts.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Retrieve structured artifacts from a stored assay
#'
#' Returns the opaque, provider-neutral artifact list persisted when a
#' structured step result was materialized. The selected assay must be stored;
#' this accessor never replays a virtual transformation or provider. Assays
#' without structured artifacts return an empty list. Unknown assays, virtual
#' targets, and stored assays with malformed non-list artifact metadata cause
#' an error.
#'
#' @param object A `ProBatchFeatures` object.
#' @param assay Stored assay identifier. `NULL` selects the current assay.
#'
#' @return The persisted artifact list, unchanged, or `list()` when the stored
#'   assay has no structured artifacts.
#' @export
#' @md
#' @examples
#' input <- matrix(
#'     1:4,
#'     nrow = 2,
#'     dimnames = list(c("f1", "f2"), c("s1", "s2"))
#' )
#' pbf <- ProBatchFeatures(input, name = "raw")
#' pb_step_artifacts(pbf)
pb_step_artifacts <- function(object, assay = NULL) {
    if (!methods::is(object, "ProBatchFeatures")) {
        stop("`object` must be a ProBatchFeatures object.", call. = FALSE)
    }
    if (is.null(assay)) {
        assay <- pb_current_assay(object)
        if (length(assay) != 1L || is.na(assay) || !nzchar(assay)) {
            stop(
                "`object` has no current stored assay.",
                call. = FALSE
            )
        }
    } else if (
        !is.character(assay) ||
            length(assay) != 1L ||
            is.na(assay) ||
            !nzchar(assay)
    ) {
        stop(
            "`assay` must be NULL or one non-missing, non-empty character value.",
            call. = FALSE
        )
    }

    if (!(assay %in% names(object))) {
        operation_log <- get_operation_log(object)
        virtual <- nrow(operation_log) &&
            any(as.character(operation_log$to) == assay)
        if (virtual) {
            stop(
                "Assay '",
                assay,
                "' is virtual and has not been materialized.",
                call. = FALSE
            )
        }
        stop("Assay '", assay, "' is not stored in the object.", call. = FALSE)
    }

    stored <- .pb_read_step_artifact_metadata(
        object[[assay]],
        paste0("Stored assay '", assay, "'")
    )
    stored$artifacts
}
