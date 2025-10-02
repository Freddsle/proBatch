.onLoad <- function(libname, pkgname) {
    safe_register <- function(name, fun) {
        try(pb_register_step(name, fun), silent = TRUE)
    }

    # existing
    safe_register("combat",   .combat_matrix_step)
    safe_register("limmaRBE", .removeBatchEffect_matrix_step)

    # NEW: only expose if the package is available
    if (.pb_requireNamespace("BERT")) {
        safe_register("BERT", .bert_matrix_step)  # canonical name
        safe_register("bert", .bert_matrix_step)  # convenience alias (lowercase)
    }
}

# internal helper to assert package availability without attaching it
.pb_requireNamespace <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf(
            "Package '%s' is required for this operation. Install via BiocManager::install('%s').",
            pkg, pkg
        ), call. = FALSE)
    }
    invisible(TRUE)
}