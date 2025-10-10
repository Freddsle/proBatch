expect_matrix_like <- function(actual, template, storage_mode = "double") {
    expect_true(is.matrix(actual))
    expect_identical(dim(actual), dim(template))
    expect_identical(dimnames(actual), dimnames(template))
    if (!is.null(storage_mode)) {
        expect_identical(storage.mode(actual), storage_mode)
    }
}

local_fake_normae_core <- function(fake_core) {
    caller_env <- parent.frame()
    testthat::local_mocked_bindings(
        .run_normae_core = fake_core,
        .normae_prepare_python = function(...) "/mock/python",
        .package = "proBatch",
        .env = caller_env
    )
}

test_that("correct_with_NormAE(wide): forwards injection order and preserves dimnames (mocked)", {
    m <- matrix(
        c(
            10, 20, 30,
            40, 50, 60
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("feat1", "feat2"), c("s1", "s2", "s3"))
    )

    sa <- data.frame(
        FullRunName = c("s3", "s1", "s2"),
        MS_batch = c("B2", "B1", "B1"),
        InjectionOrder = c(3L, 1L, 2L),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    fake_core <- function(data_matrix, sample_annotation, sample_id_col,
                          batch_col, inj_order_col, qc_col_name,
                          python_env, conda_env, normae_args, ...) {
        captured$data_matrix <- data_matrix
        captured$sample_annotation <- sample_annotation
        captured$sample_id_col <- sample_id_col
        captured$batch_col <- batch_col
        captured$inj_order_col <- inj_order_col
        captured$qc_col_name <- qc_col_name
        captured$normae_args <- normae_args
        storage.mode(data_matrix) <- "double"
        data_matrix + 5
    }

    local_fake_normae_core(fake_core)

    out <- proBatch::correct_with_NormAE(
        x = m,
        sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        inj_order_col = "InjectionOrder",
        format = "wide"
    )

    expect_matrix_like(out, m)
    expect_equal(out, m + 5, tolerance = 1e-12)
    expect_identical(captured$inj_order_col, "InjectionOrder")
    expect_null(captured$qc_col_name)
    expect_equal(captured$normae_args, list())
    expect_identical(colnames(captured$data_matrix), colnames(m))
})

test_that("correct_with_NormAE(long): adds preBatchCorr_* and respects keep_all = 'minimal' (mocked)", {
    df <- data.frame(
        peptide_group_label = rep(c("p1", "p2"), each = 3),
        FullRunName         = rep(c("s1", "s2", "s3"), times = 2),
        Intensity           = as.numeric(1:6),
        stringsAsFactors    = FALSE
    )

    sa <- data.frame(
        FullRunName = c("s3", "s1", "s2"),
        MS_batch = c("B2", "B1", "B1"),
        InjectionOrder = c(3L, 1L, 2L),
        stringsAsFactors = FALSE
    )

    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 2
    }

    local_fake_normae_core(fake_core)

    out <- correct_with_NormAE(
        x = df,
        sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity",
        batch_col = "MS_batch",
        inj_order_col = "InjectionOrder",
        format = "long",
        keep_all = "minimal"
    )

    expect_true(is.data.frame(out))
    expect_true(all(c(
        "peptide_group_label", "FullRunName",
        "Intensity", "preBatchCorr_Intensity"
    ) %in% names(out)))

    idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        interaction(out$peptide_group_label, out$FullRunName)
    )
    expect_equal(out$preBatchCorr_Intensity[idx], df$Intensity)
    expect_equal(out$Intensity[idx], df$Intensity + 2)
})

call_mocked_run_normae_core <- function(fake_py, fake_align, fake_system2, ...) {
    run_normae_core <- proBatch:::`.run_normae_core`
    patched_env <- new.env(parent = environment(run_normae_core))
    patched_env$.normae_prepare_python <- function(...) fake_py
    patched_env$.align_sample_annotation <- fake_align
    patched_env$system2 <- fake_system2
    environment(run_normae_core) <- patched_env
    run_normae_core(...)
}

test_that(".normae_matrix_step rejects matrices with NAs", {
    m <- matrix(
        c(1, 2, NA, 4),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2"),
        MS_batch = c("B1", "B2"),
        stringsAsFactors = FALSE
    )

    expect_error(
        proBatch:::`.normae_matrix_step`(
            data_matrix = m,
            sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch"
        ),
        "NormAE requires no NAs",
        fixed = TRUE
    )
})

test_that(".run_normae_core builds CLI command and returns cleaned matrix (mocked system2)", {
    data_matrix <- matrix(
        as.numeric(1:9),
        nrow = 3,
        dimnames = list(paste0("feat", 1:3), paste0("run", 1:3))
    )

    sample_annotation <- data.frame(
        FullRunName = paste0("run", c(3, 1, 2)), # shuffled order
        MS_batch = c("B2", "B1", "B1"),
        InjectionOrder = c("3", "1", "2"),
        QCFlag = c(TRUE, FALSE, TRUE),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    fake_py <- "/mock/python"

    fake_system2 <- function(command, args, stdout, stderr) {
        captured$command <- command
        captured$args <- args
        sample_csv <- args[match("--sample_csv", args) + 1L]
        meta_csv <- args[match("--meta_csv", args) + 1L]
        captured$sample_csv <- utils::read.csv(sample_csv, stringsAsFactors = FALSE, check.names = FALSE)
        captured$meta_csv <- utils::read.csv(meta_csv, row.names = 1, check.names = FALSE)
        out_dir <- args[match("--output_dir", args) + 1L]
        clean_matrix <- data_matrix + 100
        utils::write.csv(clean_matrix, file.path(out_dir, "X_clean.csv"), quote = TRUE)
        structure("NormAE CLI mock", status = 0L)
    }

    result <- call_mocked_run_normae_core(
        fake_py = fake_py,
        fake_align = function(sample_annotation, sample_ids, sample_id_col) {
            sample_annotation[match(sample_ids, sample_annotation[[sample_id_col]]), , drop = FALSE]
        },
        fake_system2 = fake_system2,
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        inj_order_col = "InjectionOrder",
        qc_col_name = "QCFlag",
        normae_args = list(
            batch_size = 32,
            enc_hiddens = c(256, 128),
            log_transform = TRUE,
            meta_csv = "ignored"
        )
    )

    expect_matrix_like(result, data_matrix)
    expect_equal(result, data_matrix + 100)
    expect_equal(captured$command, fake_py)
    expect_identical(captured$args[1:2], c("-m", "normae.cli"))
    expect_true("--meta_csv" %in% captured$args)
    expect_true("--sample_csv" %in% captured$args)
    expect_true("--output_dir" %in% captured$args)
    expect_true("--batch_indicator_col" %in% captured$args)
    expect_true("--order_indicator_col" %in% captured$args)
    expect_true("--qc_indicator_col" %in% captured$args)
    expect_true("--qc_indicator_value" %in% captured$args)

    expect_identical(
        captured$sample_csv$FullRunName,
        paste0("run", 1:3)
    )
    expect_equal(captured$sample_csv$InjectionOrder, c(1, 2, 3))
    expect_identical(captured$sample_csv$`..normae_qc`, c("Subject", "QC", "QC"))

    idx_order <- match("--order_indicator_col", captured$args)
    expect_identical(captured$args[idx_order + 1L], "InjectionOrder")

    idx_qc <- match("--qc_indicator_col", captured$args)
    expect_identical(captured$args[idx_qc + 1L], "..normae_qc")

    idx_batch <- match("--batch_size", captured$args)
    expect_false(is.na(idx_batch))
    expect_identical(captured$args[idx_batch + 1L], "32")

    idx_enc <- match("--enc_hiddens", captured$args)
    expect_false(is.na(idx_enc))
    expect_identical(captured$args[idx_enc + 1L], "256")
    expect_identical(captured$args[idx_enc + 2L], "128")

    idx_log <- match("--log_transform", captured$args)
    expect_false(is.na(idx_log))
    expect_identical(captured$args[idx_log + 1L], "True")

    expect_true("--mz_row" %in% captured$args)
    expect_true("--rt_row" %in% captured$args)
    expect_true("--early_stop" %in% captured$args)
    expect_false("ignored" %in% captured$args)
})

test_that(".run_normae_core derives default injection order when none supplied", {
    data_matrix <- matrix(
        as.numeric(1:6),
        nrow = 2,
        dimnames = list(paste0("feat", 1:2), paste0("s", 1:3))
    )

    sample_annotation <- data.frame(
        FullRunName = paste0("s", c(2, 1, 3)),
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    captured <- new.env(parent = emptyenv())
    fake_py <- "/mock/python"

    fake_system2 <- function(command, args, stdout, stderr) {
        captured$args <- args
        sample_csv <- args[match("--sample_csv", args) + 1L]
        captured$sample_csv <- utils::read.csv(sample_csv, stringsAsFactors = FALSE, check.names = FALSE)
        out_dir <- args[match("--output_dir", args) + 1L]
        clean_matrix <- data_matrix + 50
        utils::write.csv(clean_matrix, file.path(out_dir, "X_clean.csv"), quote = TRUE)
        structure("NormAE CLI mock", status = 0L)
    }

    result <- call_mocked_run_normae_core(
        fake_py = fake_py,
        fake_align = function(sample_annotation, sample_ids, sample_id_col) {
            sample_annotation[match(sample_ids, sample_annotation[[sample_id_col]]), , drop = FALSE]
        },
        fake_system2 = fake_system2,
        data_matrix = data_matrix,
        sample_annotation = sample_annotation,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        inj_order_col = NULL,
        qc_col_name = NULL,
        normae_args = list()
    )

    expect_matrix_like(result, data_matrix)
    expect_equal(result, data_matrix + 50)

    idx_order <- match("--order_indicator_col", captured$args)
    expect_identical(captured$args[idx_order + 1L], "..normae_order")
    expect_equal(captured$sample_csv$`..normae_order`, c(1, 2, 3))
})

test_that(".normae_format_cli_args skips reserved keys, injects defaults, and coerces values", {
    args <- proBatch:::`.normae_format_cli_args`(list(
        batch_size = 128,
        enc_hiddens = list(300, 200),
        log_transform = FALSE,
        meta_csv = "skip"
    ))

    expect_equal(
        args,
        c(
            "--batch_size", "128",
            "--enc_hiddens", "300", "200",
            "--log_transform", "False",
            "--mz_row", "",
            "--rt_row", "",
            "--early_stop", "False"
        )
    )
})

test_that(".normae_qc_mask handles logical, factor, and character inputs", {
    expect_identical(
        proBatch:::`.normae_qc_mask`(c(TRUE, FALSE, NA)),
        c(TRUE, FALSE, FALSE)
    )

    expect_identical(
        proBatch:::`.normae_qc_mask`(factor(c("QC", "Subject", NA))),
        c(TRUE, FALSE, FALSE)
    )

    expect_identical(
        proBatch:::`.normae_qc_mask`(c("QC", "quality_control", "random")),
        c(TRUE, TRUE, FALSE)
    )
})

test_that("correct_with_NormAE runs the NormAE CLI when available (integration)", {
    skip_on_cran()
    skip_if(nzchar(Sys.getenv("BBS_HOME")), "Skipping on Bioconductor build system")
    skip_if(
        tolower(Sys.getenv("BIOCONDUCTOR_DOCKER", "")) %in% c("true", "1"),
        "Skipping on Bioconductor docker checks"
    )
    skip_if_not_installed("reticulate")

    python_override <- NULL
    preconfigured_python <- Sys.getenv("RETICULATE_PYTHON", "")
    if (nzchar(preconfigured_python)) {
        python_bin <- preconfigured_python
    } else {
        python_bin <- Sys.which("python")
    }
    if (nzchar(python_bin)) {
        cli_probe <- tryCatch(
            system2(
                python_bin,
                c("-c", "import normae.cli"),
                stdout = TRUE,
                stderr = TRUE
            ),
            error = identity
        )
        if (!inherits(cli_probe, "error") && is.null(attr(cli_probe, "status"))) {
            python_override <- python_bin
        }
    }

    old_conda <- NULL
    if (is.null(python_override)) {
        conda_bin <- Sys.which("mamba")
        if (!nzchar(conda_bin)) {
            conda_bin <- Sys.which("conda")
        }
        skip_if(!nzchar(conda_bin), "No mamba/conda binary available for NormAE env creation")
        old_conda <- getOption("reticulate.conda_binary", NULL)
        options(reticulate.conda_binary = conda_bin)
        on.exit(options(reticulate.conda_binary = old_conda), add = TRUE)
    }

    set.seed(123)
    m <- matrix(
        runif(12, min = 100, max = 1000),
        nrow = 3,
        dimnames = list(
            paste0("feat", 1:3),
            paste0("sample", 1:4)
        )
    )

    sa <- data.frame(
        FullRunName = paste0("sample", 1:4),
        MS_batch = rep(c("B1", "B2"), each = 2),
        InjectionOrder = 1:4,
        stringsAsFactors = FALSE
    )

    env_error_patterns <- paste(
        c(
            "Python module 'normae'",
            "Python module 'normae.cli'",
            "No module named 'normae'",
            "conda binary",
            "CondaEnvException",
            "Miniconda",
            "conda environment 'normae' not found"
        ),
        collapse = "|"
    )

    result <- tryCatch(
        correct_with_NormAE(
            x = m,
            sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            inj_order_col = "InjectionOrder",
            format = "wide",
            python_env = python_override,
            normae_args = list()
        ),
        error = function(e) {
            msg <- conditionMessage(e)
            if (grepl(env_error_patterns, msg, ignore.case = TRUE)) {
                if (!is.null(python_override)) {
                    skip(paste("NormAE CLI unavailable in system python:", msg))
                }
                skip(paste("NormAE environment not ready:", msg))
            }
            stop(msg, call. = FALSE)
        }
    )

    expect_matrix_like(result, m)
    expect_false(isTRUE(all.equal(result, m)))
    expect_true(all(is.finite(result)))
})
