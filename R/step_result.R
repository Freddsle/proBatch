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
    assay_metadata$pb_step_artifacts <- artifacts
    S4Vectors::metadata(assay) <- assay_metadata
    assay
}
