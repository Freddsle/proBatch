# -------------------------
# Optional integration test (only if BERT is available)
# -------------------------
test_that("correct_with_BERT(wide): preserves dimnames and returns numeric matrix (mocked)", {
    testthat::skip_if_not_installed("BERT")
    m <- matrix(
        c(
            1, NA, 3,
            4, 5, 6
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )

    sa <- data.frame(
        FullRunName = c("s2", "s1", "s3"), # shuffled to test alignment
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    # Fake core that adds +10 to all numeric cells and returns same shape
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 10
    }

    local_mocked_bindings(
        .pb_requireNamespace = function(pkg) invisible(TRUE),
        .bert_matrix_step    = fake_core,
        .package             = "proBatch"
    )

    out <- proBatch::correct_with_BERT(
        x = m, sample_annotation = sa,
        sample_id_col = "FullRunName",
        batch_col = "MS_batch",
        format = "wide"
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
    expect_identical(storage.mode(out), "double")
    expect_equal(out, m + 10, tolerance = 1e-12)
})

test_that("correct_with_BERT(long): adds preBatchCorr_* and respects keep_all = 'minimal' (mocked)", {
    testthat::skip_if_not_installed("BERT")
    # Toy long with two features × three samples (complete)
    df <- data.frame(
        peptide_group_label = rep(c("p1", "p2"), each = 3),
        FullRunName         = rep(c("s1", "s2", "s3"), times = 2),
        Intensity           = as.numeric(1:6),
        stringsAsFactors    = FALSE
    )

    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3"),
        MS_batch = c("B1", "B1", "B2"),
        stringsAsFactors = FALSE
    )

    # Core returns +1 everywhere
    fake_core <- function(data_matrix, ...) {
        storage.mode(data_matrix) <- "double"
        data_matrix + 1
    }

    local_mocked_bindings(
        .pb_requireNamespace = function(pkg) invisible(TRUE),
        .bert_matrix_step    = fake_core,
        .package             = "proBatch"
    )

    out <- correct_with_BERT(
        x = df, sample_annotation = sa,
        feature_id_col = "peptide_group_label",
        sample_id_col = "FullRunName",
        measure_col = "Intensity",
        batch_col = "MS_batch",
        format = "long",
        keep_all = "minimal"
    )
    expect_true(is.data.frame(out))

    # Minimal set of columns (order not enforced)
    expect_true(all(c(
        "peptide_group_label", "FullRunName", "Intensity",
        "preBatchCorr_Intensity"
    ) %in% names(out)))

    # preBatchCorr_* must match the original measure by (feature, sample)
    idx <- match(
        interaction(df$peptide_group_label, df$FullRunName),
        interaction(out$peptide_group_label, out$FullRunName)
    )
    expect_equal(out$preBatchCorr_Intensity[idx], df$Intensity)

    # Corrected measure is +1 vs original
    expect_equal(out$Intensity[idx], df$Intensity + 1)
})

test_that(".bert_matrix_step: rejects non-numeric and <2 features (mocked/no BERT)", {
    testthat::skip_if_not_installed("BERT")
    m_bad <- matrix(
        c("a", "b", "c", "d"),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
    sa <- data.frame(FullRunName = c("s1", "s2"), MS_batch = c("B1", "B2"))

    expect_error(
        proBatch:::`.bert_matrix_step`(
            data_matrix = m_bad, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch"
        ),
        "Input must be coercible to a numeric matrix for BERT correction",
        fixed = TRUE
    )

    m_tiny <- matrix(1:2, nrow = 1, dimnames = list("f1", c("s1", "s2")))
    expect_error(
        proBatch:::`.bert_matrix_step`(
            data_matrix = m_tiny, sample_annotation = sa,
            sample_id_col = "FullRunName", batch_col = "MS_batch"
        ),
        " at least two not-NA features.",
        ignore.case = TRUE
    )
})

test_that(".run_BERT_core: fails if batch column missing; aligns SA to matrix columns", {
    testthat::skip_if_not_installed("BERT")
    m <- matrix(
        c(
            1.0, 2.0, 3.0, 4.0,
            4.0, 5.0, 6.0, 7.0
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )
    # Missing batch column -> error
    sa_bad <- data.frame(FullRunName = c("s3", "s2", "s1", "s4"))

    expect_error(
        proBatch:::`.run_BERT_core`(
            data_matrix = m, sample_annotation = sa_bad,
            sample_id_col = "FullRunName", batch_col = "MS_batch",
            bert_method = "limma"
        ),
        "Batch column is not present",
        fixed = TRUE
    )

    # Alignment check with shuffled SA
    sa <- data.frame(
        FullRunName = c("s3", "s2", "s1", "s4"), # shuffled on purpose
        MS_batch = c(1, 2, 1, 2),
        stringsAsFactors = FALSE
    )
    out <- proBatch:::`.run_BERT_core`(
        data_matrix = m, sample_annotation = sa,
        sample_id_col = "FullRunName", batch_col = "MS_batch",
        bert_method = "limma"
    )
    expect_true(is.matrix(out))
    expect_identical(rownames(out), rownames(m))
    expect_identical(colnames(out), colnames(m))
})

test_that("correct_with_BERT(long/wide): format argument checking and plumbing (mocked)", {
    testthat::skip_if_not_installed("BERT")
    sa <- data.frame(FullRunName = c("s1", "s2"), MS_batch = c("B1", "B2"))

    # format='wide' but x not a matrix -> error
    expect_error(
        proBatch::correct_with_BERT(
            x = data.frame(x = 1:2), sample_annotation = sa,
            format = "wide"
        ),
        "requires a numeric matrix",
        ignore.case = TRUE
    )

    # format='long' but x not a data.frame -> error
    expect_error(
        proBatch::correct_with_BERT(
            x = matrix(1, 1, 1), sample_annotation = sa,
            format = "long"
        ),
        "requires a data.frame",
        ignore.case = TRUE
    )
})

test_that("Integration: .bert_matrix_step runs with BERT on tiny data (method='limma')", {
    testthat::skip_if_not_installed("BERT")

    m <- matrix(
        c(
            1, NA, 3, 4,
            4, 5, 6, 7
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    )
    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3", "s4"),
        MS_batch = c("B1", "B1", "B2", "B2"),
        Condition = c("A", "A", "B", "B"), # covariate
        stringsAsFactors = FALSE
    )

    out <- proBatch:::`.bert_matrix_step`(
        data_matrix       = m,
        sample_annotation = sa,
        sample_id_col     = "FullRunName",
        batch_col         = "MS_batch",
        covariates_cols   = "Condition",
        bert_method       = "limma", # avoids ComBat priors & works with NA per BERT
        combatmode        = 1
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
    expect_equal(storage.mode(out), "double")
})


# -------------------------
# Optional integration test (only if PLSDAbatch is available)
# -------------------------
test_that("correct_with_PLSDAbatch(wide)", {
    testthat::skip_if_not_installed("PLSDAbatch")
    m <- matrix(
        c(
            1, 2, 3, 4, 1,
            4, 5, 6, 7, 4
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(
            c("f1", "f2"),
            c("s1", "s2", "s3", "s4", "s5")
        )
    )

    sa <- data.frame(
        FullRunName = c("s1", "s2", "s3", "s4", "s5"),
        MS_batch = c("B1", "B1", "B2", "B2", "B1"),
        Condition = c("A", "A", "B", "B", "A"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        out <- correct_with_PLSDA_batch(
            x = m, sample_annotation = sa,
            sample_id_col = "FullRunName",
            batch_col = "MS_batch",
            format = "wide",
            effect_col = "Condition",
        ),
        "Some treatment levels have fewer than 3 samples"
    )

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(m))
    expect_identical(dimnames(out), dimnames(m))
    expect_identical(storage.mode(out), "double")
})

test_that("correct_with_PLSDA_batch(wide, mocked): returns double matrix and preserves dimnames", {

  testthat::skip_if_not_installed("PLSDAbatch")
  # Small toy matrix (features x samples)
  m <- matrix(
    c(
      1, 2, 3, 4, 5,
      4, 5, 6, 7, 8
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("s", 1:5))
  )

  # Unbalanced design (B2 has only 'B'): forces balance=FALSE in "auto"
  sa <- data.frame(
    FullRunName = colnames(m),
    MS_batch    = c("B1", "B1", "B1", "B2", "B2"),
    Condition   = c("A", "A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  # Capture arguments reaching PLSDAbatch::PLSDA_batch
  capt <- new.env(parent = emptyenv())

  fake_plsda_batch <- function(...) {
    args <- list(...)
    X <- args$X # samples x features
    # record what the wrapper passed
    capt$last_balance        <- isTRUE(args$balance)
    capt$last_keepX_trt      <- args$keepX.trt
    capt$last_Y_trt_is_null  <- is.null(args$Y.trt)

    nbat <- if (!is.null(args$ncomp.bat)) args$ncomp.bat else 2L
    # minimal structure used by .select_ncomp_bat()
    list(
      X.nobatch = X, # no-op correction; we only test shape/invariants
      explained_variance.bat = list(
        X = rep(0, nbat),
        Y = rep(1 / nbat, nbat)
      )
    )
  }

  fake_plsda <- function(X, Y, ncomp) {
    list(prop_expl_var = list(
      X = rep(0, ncomp),
      Y = rep(1 / ncomp, ncomp)
    ))
  }

  fake_tune <- function(X, Y, ncomp, test.keepX, ...) {
    list(choice.keepX = rep(11L, ncomp))
  }

  # Mock PLSDAbatch & mixOmics entrypoints
  testthat::local_mocked_bindings(
    PLSDA_batch = fake_plsda_batch,
    .env = asNamespace("PLSDAbatch")
  )
  testthat::local_mocked_bindings(
    plsda = fake_plsda,
    `tune.splsda` = fake_tune,
    .env = asNamespace("mixOmics")
  )

  # 1) basic wide run (auto-select ncomp_trt/bat; unbalanced=FALSE expected)
  out <- correct_with_PLSDA_batch(
    x = m, sample_annotation = sa,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    effect_col = "Condition",
    format = "wide"
  )

  expect_true(is.matrix(out))
  expect_identical(dim(out), dim(m))
  expect_identical(dimnames(out), dimnames(m))
  expect_identical(storage.mode(out), "double")
  # auto-balance must flip to FALSE for zero cells in table(batch, trt)
  expect_identical(capt$last_balance, FALSE)

  # 2) effect_col = NULL (should pass Y.trt = NULL into PLSDAbatch)
  out2 <- correct_with_PLSDA_batch(
    x = m, sample_annotation = sa,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    effect_col = NULL,
    format = "wide"
  )
  expect_true(is.matrix(out2))
  expect_true(isTRUE(capt$last_Y_trt_is_null))

  # 3) sPLSDA path with tuning (keepX_trt inferred via tune.splsda)
  out3 <- correct_with_PLSDA_batch(
    x = m, sample_annotation = sa,
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    effect_col = "Condition",
    ncomp_trt = 2L,
    run_splsda = TRUE,
    format = "wide"
  )
  expect_true(is.matrix(out3))
  expect_identical(capt$last_keepX_trt, c(11L, 11L))
})

test_that("correct_with_PLSDA_batch(wide): basic validation errors are informative", {

  testthat::skip_if_not_installed("PLSDAbatch")
  m <- matrix(c(1, 2, NA, 4), nrow = 2,
              dimnames = list(c("f1", "f2"), c("s1", "s2")))
  sa <- data.frame(
    FullRunName = colnames(m),
    MS_batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )

  # Mock PLSDAbatch to avoid loading the package; we won't reach it anyway
  testthat::local_mocked_bindings(
    PLSDA_batch = function(...) stop("should not be called"),
    .env = asNamespace("PLSDAbatch")
  )
  testthat::local_mocked_bindings(
    plsda = function(...) stop("should not be called"),
    `tune.splsda` = function(...) stop("should not be called"),
    .env = asNamespace("mixOmics")
  )

  expect_error(
    correct_with_PLSDA_batch(
      x = m, sample_annotation = sa,
      sample_id_col = "FullRunName",
      batch_col = "MS_batch",
      effect_col = "Condition",
      format = "wide"
    ),
    "requires no NAs", ignore.case = TRUE
  )

  expect_error(
    correct_with_PLSDA_batch(
      x = as.data.frame(m), sample_annotation = sa,
      sample_id_col = "FullRunName",
      batch_col = "MS_batch",
      effect_col = "Condition",
      format = "wide"
    ),
    "requires a numeric matrix", ignore.case = TRUE
  )

  expect_error(
    correct_with_PLSDA_batch(
      x = matrix(c(1, 2, 3, 4), 2,
                 dimnames = list(c("f1", "f2"), c("s1", "s2"))),
      sample_annotation = sa[, "FullRunName", drop = FALSE], # drop batch col
      sample_id_col = "FullRunName",
      batch_col = "MS_batch",
      format = "wide"
    ),
    "Batch column is not present", fixed = TRUE
  )
})

test_that("correct_with_PLSDA_batch(long, mocked): round-trip keeps rows and adds preBatchCorr_*", {

  testthat::skip_if_not_installed("PLSDAbatch")
  # make a 2x3 matrix and longify it
  m <- matrix(
    c(
      1, 2, 3,
      4, 5, 6
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
  )
  df_long <- data.frame(
    peptide_group_label = rep(rownames(m), times = ncol(m)),
    FullRunName = rep(colnames(m), each = nrow(m)),
    Intensity = as.numeric(m),
    extra_col = "keepme",
    stringsAsFactors = FALSE
  )
  sa <- data.frame(
    FullRunName = colnames(m),
    MS_batch = c("B1", "B1", "B2"),
    Condition = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  # very light mocks (no-ops)
  testthat::local_mocked_bindings(
    PLSDA_batch = function(X, Y.trt, Y.bat, ncomp.trt, ncomp.bat, ..., balance = TRUE) {
      nb <- if (is.null(ncomp.bat)) 1L else ncomp.bat
      list(
        X.nobatch = X,
        explained_variance.bat = list(X = rep(0, nb), Y = rep(1/nb, nb))
      )
    },
    .env = asNamespace("PLSDAbatch")
  )
  testthat::local_mocked_bindings(
    plsda = function(X, Y, ncomp) list(prop_expl_var = list(X = rep(0, ncomp), Y = rep(1/ncomp, ncomp))),
    `tune.splsda` = function(...) list(choice.keepX = 5L),
    .env = asNamespace("mixOmics")
  )

  out <- correct_with_PLSDA_batch(
    x = df_long,
    sample_annotation = sa,
    sample_id_col = "FullRunName",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    batch_col = "MS_batch",
    effect_col = "Condition",
    format = "long"
  )

  expect_true(is.data.frame(out))
  expect_true(all(c("peptide_group_label", "FullRunName", "Intensity") %in% names(out)))
  precol <- "preBatchCorr_Intensity"
  expect_true(precol %in% names(out))
  # row count is preserved (features * samples that were present)
  expect_equal(nrow(out), nrow(df_long))

  # minimal columns request
  out_min <- correct_with_PLSDA_batch(
    x = df_long,
    sample_annotation = sa,
    sample_id_col = "FullRunName",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    batch_col = "MS_batch",
    effect_col = "Condition",
    format = "long",
    keep_all = "minimal"
  )
  expect_true(setequal(
    names(out_min),
    c("FullRunName", "peptide_group_label", "Intensity", precol)
  ))
})


test_that("correct_with_PLSDA_batch(wide, integration) runs with PLSDAbatch installed", {
  testthat::skip_if_not_installed("PLSDAbatch")

  m <- matrix(
    c(
      1, 2, 3, 4, 1,
      4, 5, 6, 7, 4
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("s", 1:5))
  )
  sa <- data.frame(
    FullRunName = colnames(m),
    MS_batch = c("B1", "B1", "B2", "B2", "B1"),
    Condition = c("A", "A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    out <- correct_with_PLSDA_batch(
      x = m, sample_annotation = sa,
      sample_id_col = "FullRunName",
      batch_col = "MS_batch",
      effect_col = "Condition",
      format = "wide"
    ),
    "fewer than 3 samples", ignore.case = TRUE
  )
  expect_true(is.matrix(out))
  expect_identical(dim(out), dim(m))
  expect_identical(dimnames(out), dimnames(m))
  expect_identical(storage.mode(out), "double")
})

test_that("correct_with_PLSDA_batch(long, integration) runs with PLSDAbatch installed", {
  testthat::skip_if_not_installed("PLSDAbatch")

  m <- matrix(
    c(1,2,3, 4,5,6),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("f1","f2"), c("s1","s2","s3"))
  )
  df_long <- data.frame(
    peptide_group_label = rep(rownames(m), times = ncol(m)),
    FullRunName = rep(colnames(m), each = nrow(m)),
    Intensity = as.numeric(m),
    stringsAsFactors = FALSE
  )
  sa <- data.frame(
    FullRunName = colnames(m),
    MS_batch = c("B1", "B1", "B2"),
    Condition = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  out <- correct_with_PLSDA_batch(
    x = df_long, sample_annotation = sa,
    sample_id_col = "FullRunName",
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    batch_col = "MS_batch",
    effect_col = "Condition",
    format = "long"
  )
  expect_true(is.data.frame(out))
  expect_true(all(c("peptide_group_label","FullRunName","Intensity","preBatchCorr_Intensity") %in% names(out)))
})

test_that("ProBatchFeatures smoke-test (skipped if class/methods unavailable)", {
  testthat::skip_if_not_installed("PLSDAbatch")

  m <- matrix(
    c(1,2,3, 4,5,6),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("f1","f2"), c("s1","s2","s3"))
  )
  sa <- S4Vectors::DataFrame(
    FullRunName = colnames(m),
    MS_batch = c("B1","B1","B2"),
    Condition = c("A","A","B")
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = m),
    colData = sa
  )
  qf <- QFeatures::QFeatures(list(`protein::raw` = se))

  pb <- try(methods::as(qf, "ProBatchFeatures"), silent = TRUE)
  if (inherits(pb, "try-error")) testthat::skip("Coercion to ProBatchFeatures not available in this branch")

  # Use the matrix from the assay to exercise the (wide) path with SA from pb
  out <- correct_with_PLSDA_batch(
    x = SummarizedExperiment::assay(QFeatures::assay(pb, 1L)),
    sample_annotation = as.data.frame(SummarizedExperiment::colData(QFeatures::assay(pb, 1L))),
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    effect_col = "Condition",
    format = "wide"
  )
  expect_true(is.matrix(out))
})