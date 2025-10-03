.onLoad <- function(libname, pkgname) {
    safe_register <- function(name, fun) {
        try(pb_register_step(name, fun), silent = TRUE)
    }

    # existing
    safe_register("combat", .combat_matrix_step)
    safe_register("limmaRBE", .removeBatchEffect_matrix_step)
    # mComBat uses the same func, but with use_mComBat = TRUE
    safe_register("mComBat", .mComBat_matrix_step)

    # Only expose if the package is available
    if (.pb_requireNamespace("BERT", only_info = TRUE)) {
        safe_register("BERT", .bert_matrix_step) # canonical name
        safe_register("bert", .bert_matrix_step) # convenience alias (lowercase)
    }
}

# internal helper to assert package availability without attaching it
.pb_requireNamespace <- function(pkg, only_info = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (only_info) {
            packageStartupMessage(
                sprintf(
                    "Package '%s' is not installed. If this functionality is required, please install via BiocManager::install('%s').",
                    pkg, pkg
                )
            )
            return(invisible(FALSE))
        }
        stop(sprintf(
            "Package '%s' is required for this operation. Install via BiocManager::install('%s').",
            pkg, pkg
        ), call. = FALSE)
    }
    invisible(TRUE)
}
