data(example_proteome, package = "proBatch")
data(example_sample_annotation, package = "proBatch")

# Default arguments for calculate_feature_CV
default_args <- list(
    df_long            = example_proteome,
    sample_annotation  = example_sample_annotation,
    sample_id_col      = NULL, # let calculate_feature_CV pick defaults
    feature_id_col     = "peptide_group_label",
    measure_col        = "Intensity",
    batch_col          = NULL,
    biospecimen_id_col = "EarTag",
    unlog              = FALSE
)

# Helper to call with overrides
calc_cv <- function(...) {
    args <- modifyList(default_args, list(...))
    do.call(calculate_feature_CV, args)
}

test_that("basic output structure without batch, default args", {
    cv_df <- calc_cv()
    expect_equal(names(cv_df), c("peptide_group_label", "CV_total"))
    expect_type(cv_df$CV_total, "double")
    expect_true(all(cv_df$CV_total >= 0))
})

test_that("error on missing biospecimen_id_col", {
    expect_error(
        calc_cv(biospecimen_id_col = "NON_EXISTENT"),
        "biospecimen ID, indicating replicates, is not in the data"
    )
})

test_that("warning if no biospecimen_id_col provided", {
    expect_warning(
        cv_df <- calc_cv(biospecimen_id_col = NULL),
        "considering all samples as replicates!"
    )
    expect_true("CV_total" %in% names(cv_df))
})

test_that("filter out features with ≤2 measurements", {
    df2 <- head(example_proteome, 100)
    cv2 <- calc_cv(df_long = df2)
    low_n <- df2 %>%
        dplyr::count(peptide_group_label) %>%
        dplyr::filter(n <= 2) %>%
        dplyr::pull(peptide_group_label)
    expect_false(any(cv2$peptide_group_label %in% low_n))
})

test_that("per-batch CV appears when batch_col is supplied", {
    cvb <- calc_cv(batch_col = "MS_batch")
    expect_true(all(c("CV_perBatch", "CV_total") %in% names(cvb)))
})

test_that("unlog = FALSE keeps log-scale, no unlog warning", {
    expect_silent(
        cvn <- calc_cv(unlog = FALSE)
    )
})

test_that("features with <=2 replicates are removed with warning", {
    df <- tibble::tibble(
        peptide_group_label = c(rep("pep1", 2), rep("pep2", 3)),
        Intensity           = c(1, 2, 10, 12, 11),
        FullRunName         = paste0("sample", 1:5),
        MS_batch            = rep("A", 5),
        biospec             = rep("bio1", 5)
    )
    expect_warning(
        res <- calculate_feature_CV(
            df_long            = df,
            sample_annotation  = NULL,
            sample_id_col      = "FullRunName",
            feature_id_col     = "peptide_group_label",
            measure_col        = "Intensity",
            batch_col          = "MS_batch",
            biospecimen_id_col = "biospec",
            unlog              = FALSE
        ),
        "removing those peptides"
    )
    expect_false("pep1" %in% res$peptide_group_label)
})

test_that("unlogging reproduces original CV", {
    df <- tibble::tibble(
        peptide_group_label = rep("pep", 3),
        Intensity           = c(4, 5, 6),
        FullRunName         = paste0("s", 1:3),
        biospec             = rep("bio", 3)
    )
    logged <- dplyr::mutate(df, Intensity = log2(Intensity))
    expect_warning(
        cv_logged <- calculate_feature_CV(
            df_long            = logged,
            sample_annotation  = NULL,
            sample_id_col      = "FullRunName",
            feature_id_col     = "peptide_group_label",
            measure_col        = "Intensity",
            biospecimen_id_col = "biospec",
            unlog              = TRUE
        ),
        "reversing log-transformation"
    )
    cv_raw <- calculate_feature_CV(
        df_long            = df,
        sample_annotation  = NULL,
        sample_id_col      = "FullRunName",
        feature_id_col     = "peptide_group_label",
        measure_col        = "Intensity",
        biospecimen_id_col = "biospec",
        unlog              = FALSE
    )
    expect_equal(cv_logged$CV_total, cv_raw$CV_total)
})

test_that("Step column preserved and plotting works", {
    df <- tibble::tibble(
        peptide_group_label = rep("pep", 6),
        Intensity           = c(1, 2, 1.5, 10, 11, 9),
        FullRunName         = paste0("s", 1:6),
        MS_batch            = rep(c("B1", "B2"), each = 3),
        biospec             = rep(c("bio1", "bio2"), each = 3),
        Step                = rep(c("raw", "norm"), each = 3)
    )
    cv <- calculate_feature_CV(
        df_long            = df,
        sample_annotation  = NULL,
        sample_id_col      = "FullRunName",
        feature_id_col     = "peptide_group_label",
        measure_col        = "Intensity",
        batch_col          = "MS_batch",
        biospecimen_id_col = "biospec",
        unlog              = FALSE
    )
    expect_true("Step" %in% names(cv))
    gg <- plot_CV_distr.df(cv)
    expect_s3_class(gg, "ggplot")
    expect_equal(gg$mapping$x, rlang::sym("Step"))
})

test_that("per-batch CV matches manual calculation", {
    df <- tibble::tibble(
        peptide_group_label = rep("pep", 4),
        Intensity           = c(1, 2, 5, 6),
        FullRunName         = paste0("s", 1:4),
        MS_batch            = rep(c("b1", "b2"), each = 2),
        biospec             = rep(c("bio1", "bio2"), each = 2)
    )
    cv <- calculate_feature_CV(
        df_long            = df,
        sample_annotation  = NULL,
        sample_id_col      = "FullRunName",
        feature_id_col     = "peptide_group_label",
        measure_col        = "Intensity",
        batch_col          = "MS_batch",
        biospecimen_id_col = "biospec",
        unlog              = FALSE
    )
    manual <- df %>%
        dplyr::group_by(MS_batch) %>%
        dplyr::summarise(cv = sd(Intensity) / mean(Intensity))
    expect_equal(cv$CV_perBatch[1], manual$cv[1])
    expect_equal(cv$CV_perBatch[2], manual$cv[2])
})

test_that("missing biospecimen column triggers warning", {
    expect_warning(
        res <- calculate_feature_CV(
            df_long = tibble::tibble(
                peptide_group_label = rep("pep", 3),
                Intensity           = c(1, 2, 3),
                FullRunName         = paste0("s", 1:3)
            ),
            sample_annotation = NULL,
            sample_id_col = "FullRunName",
            feature_id_col = "peptide_group_label",
            measure_col = "Intensity",
            batch_col = NULL,
            biospecimen_id_col = NULL,
            unlog = FALSE
        ),
        "considering all samples as replicates"
    )
    expect_true(all(res$CV_total >= 0))
})


# # Prepare minimal CV_df
# cv_df_min <- data.frame(
#   peptide_group_label = letters[1:10],
#   CV_total            = runif(10, 0.1, 0.5),
#   Step                = rep(c("raw","norm"), each = 5),
#   stringsAsFactors    = FALSE
# )

# test_that("returns a ggplot object", {
#   p <- plot_CV_distr.df(cv_df_min, plot_title = "Test", log_y_scale = FALSE)
#   expect_s3_class(p, "ggplot")
# })

# test_that("boxplot maps Step → x when Step present", {
#   p <- plot_CV_distr.df(cv_df_min, log_y_scale = FALSE)
#   # check mapping
#   aes_map <- layer_data(p, 1)
#   expect_true("xmin" %in% names(aes_map))  # presence of x grouping
# })

# test_that("applies log scale when requested", {
#   p_log <- plot_CV_distr.df(cv_df_min, log_y_scale = TRUE)
#   # scale_y_log10 layer must be present
#   scales <- sapply(p_log$scales$scales, class)
#   expect_true(any(grepl("ScaleContinuousPosition", scales)))
# })

# test_that("classic theme is applied when theme='classic'", {
#   p_cl <- plot_CV_distr.df(cv_df_min, theme = "classic", log_y_scale = FALSE)
#   # theme_classic adds a panel.background element of class element_blank
#   expect_true(inherits(p_cl$theme$panel.background, "element_blank"))
# })


# test_that("full pipeline returns ggplot", {
#   ggp <- plot_CV_distr(
#     df_long           = example_proteome,
#     sample_annotation = example_sample_annotation,
#     measure_col       = "Intensity",
#     batch_col         = "MS_batch",
#     biospecimen_id_col= "EarTag",
#     unlog             = TRUE,
#     plot_title        = "Full CV Test"
#   )
#   expect_s3_class(ggp, "ggplot")
#   # title must match
#   expect_equal(ggp$labels$title, "Full CV Test")
# })

# test_that("filename argument saves a file", {
#   tmpfile <- tempfile(fileext = ".png")
#   ggs <- plot_CV_distr.df(
#     CV_df     = calculate_feature_CV(
#       df_long            = example_proteome,
#       sample_annotation  = example_sample_annotation,
#       batch_col          = "MS_batch",
#       biospecimen_id_col = "EarTag"
#     ),
#     filename  = tmpfile,
#     log_y_scale = FALSE
#   )
#   expect_true(file.exists(tmpfile))
#   unlink(tmpfile)
# })
