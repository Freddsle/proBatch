.onLoad <- function(libname, pkgname) {
    register_core_step <- function(name, fun, aliases = character(),
        label = name) {
        pb_register_step(
            name = name,
            fun = fun,
            package = pkgname,
            kind = "transform",
            aliases = aliases,
            label = label
        )
    }

    register_core_step(
        "log2",
        .pb_builtin_log2_step,
        label = "Log2 transformation"
    )
    register_core_step(
        "log",
        .pb_builtin_log_step,
        label = "Log transformation"
    )
    register_core_step(
        "medianNorm",
        .pb_builtin_median_norm_step,
        label = "Median normalization"
    )
    register_core_step("combat", .combat_matrix_step)
    register_core_step("limmaRBE", .removeBatchEffect_matrix_step)
    register_core_step(
        "loessLimmaRBE",
        .loess_limmaRBE_matrix_step,
        aliases = "loesslimmarbe"
    )
    invisible()
}
