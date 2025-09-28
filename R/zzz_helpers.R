.onLoad <- function(libname, pkgname) {
    safe_register <- function(name, fun) {
        try(pb_register_step(name, fun), silent = TRUE)
    }

    safe_register("combat", .combat_matrix_step)
    safe_register("limmaRBE", .removeBatchEffect_matrix_step)
}
