.onLoad <- function(libname, pkgname) {
    safe_register <- function(name, fun) {
        try(pb_register_step(name, fun), silent = TRUE)
    }

    # existing
    safe_register("combat", .combat_matrix_step)
    safe_register("limmaRBE", .removeBatchEffect_matrix_step)
    # mComBat uses the same func, but with use_mComBat = TRUE
    safe_register("mComBat", .mComBat_matrix_step)

    # NormAE depends on Python via reticulate; only expose when infrastructure is present
    if (.pb_requireNamespace("reticulate", only_info = TRUE)) {
        safe_register("NormAE", .normae_matrix_step)
        safe_register("normae", .normae_matrix_step)
    }

    # Only expose if the package is available
    if (.pb_requireNamespace("BERT", only_info = TRUE)) {
        safe_register("BERT", .bert_matrix_step) # canonical name
        safe_register("bert", .bert_matrix_step) # convenience alias (lowercase)
    }
    if (.pb_requireNamespace("PLSDAbatch", only_info = TRUE)) {
        safe_register("PLSDAbatch", .plsda_matrix_step) # canonical name
        safe_register("plsdabatch", .plsda_matrix_step) # convenience alias (lowercase)
        safe_register("sPLSDAbatch", .splsda_matrix_step) # alias with special char
        safe_register("splsdabatch", .splsda_matrix_step) # convenience alias (lowercase)
    }

    if (.pb_requireNamespace("RUVIIIC", only_info = TRUE)) {
        safe_register("RUVIIIC", .ruviiic_matrix_step)
        safe_register("ruviiic", .ruviiic_matrix_step)
        safe_register("RUVIII_C", .ruviiic_matrix_step)
    }
}

# internal helper to assert package availability without attaching it
.pb_requireNamespace <- function(pkg, only_info = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg %in% c("reticulate")) {
            how_to_install <- "install.packages('reticulate')"
        } else if (pkg %in% c("RUVIIIC")) {
            how_to_install <- "devtools::install_github('CMRI-ProCan/RUV-III-C').\nSee https://github.com/CMRI-ProCan/RUV-III-C for more information"
        } else {
            how_to_install <- sprintf("BiocManager::install('%s')", pkg)
        }

        if (only_info) {
            packageStartupMessage(
                sprintf(
                    "Package '%s' is not installed. If this functionality is required, please install via %s.",
                    pkg, how_to_install
                )
            )
            return(invisible(FALSE))
        }
        stop(sprintf(
            "Package '%s' is required for this operation. Install via %s and try again.",
            pkg, how_to_install
        ), call. = FALSE)
    }
    invisible(TRUE)
}
