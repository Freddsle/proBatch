#' @title DEPRECATED: center_feature_batch_medians_dm
#' @description Use [center_feature_batch()] with `format="wide", stat="medians"`.
#' @inheritParams correct_batch_effects
#' @export
#' @keywords internal
center_feature_batch_medians_dm <- function(data_matrix, sample_annotation,
                                            sample_id_col = "FullRunName",
                                            batch_col = "MS_batch",
                                            feature_id_col = "peptide_group_label",
                                            measure_col = "Intensity") {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = data_matrix, sample_annotation = sample_annotation,
        format = "wide", stat = "medians",
        sample_id_col = sample_id_col, batch_col = batch_col,
        feature_id_col = feature_id_col, measure_col = measure_col
    )
}

#' @title DEPRECATED: center_feature_batch_means_dm
#' @description Use [center_feature_batch()] with `format="wide", stat="means"`.
#' @inheritParams correct_batch_effects
#' @export
#' @keywords internal
center_feature_batch_means_dm <- function(data_matrix, sample_annotation,
                                          sample_id_col = "FullRunName",
                                          batch_col = "MS_batch",
                                          feature_id_col = "peptide_group_label",
                                          measure_col = "Intensity") {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = data_matrix, sample_annotation = sample_annotation,
        format = "wide", stat = "means",
        sample_id_col = sample_id_col, batch_col = batch_col,
        feature_id_col = feature_id_col, measure_col = measure_col
    )
}

#' @title DEPRECATED: center_feature_batch_medians_df
#' @description Use [center_feature_batch()] with `format="long", stat="medians"`.
#' @inheritParams correct_batch_effects
#' @export
#' @keywords internal
center_feature_batch_medians_df <- function(df_long, sample_annotation = NULL,
                                            sample_id_col = "FullRunName",
                                            batch_col = "MS_batch",
                                            feature_id_col = "peptide_group_label",
                                            measure_col = "Intensity",
                                            keep_all = "default",
                                            no_fit_imputed = TRUE,
                                            qual_col = NULL,
                                            qual_value = NULL) {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = df_long, sample_annotation = sample_annotation,
        format = "long", stat = "medians",
        sample_id_col = sample_id_col, batch_col = batch_col,
        feature_id_col = feature_id_col, measure_col = measure_col,
        keep_all = keep_all, no_fit_imputed = no_fit_imputed,
        qual_col = qual_col, qual_value = qual_value
    )
}

#' @title DEPRECATED: center_feature_batch_means_df
#' @description Use [center_feature_batch()] with `format="long", stat="means"`.
#' @inheritParams correct_batch_effects
#' @export
#' @keywords internal
center_feature_batch_means_df <- function(df_long, sample_annotation = NULL,
                                          sample_id_col = "FullRunName",
                                          batch_col = "MS_batch",
                                          feature_id_col = "peptide_group_label",
                                          measure_col = "Intensity",
                                          keep_all = "default",
                                          no_fit_imputed = TRUE,
                                          qual_col = NULL,
                                          qual_value = NULL) {
    .Deprecated("center_feature_batch")
    center_feature_batch(
        x = df_long, sample_annotation = sample_annotation,
        format = "long", stat = "means",
        sample_id_col = sample_id_col, batch_col = batch_col,
        feature_id_col = feature_id_col, measure_col = measure_col,
        keep_all = keep_all, no_fit_imputed = no_fit_imputed,
        qual_col = qual_col, qual_value = qual_value
    )
}


#' @export
#' @rdname correct_batch_effects
correct_with_ComBat_df <- function(df_long, sample_annotation = NULL,
                                   feature_id_col = "peptide_group_label",
                                   measure_col = "Intensity",
                                   sample_id_col = "FullRunName",
                                   batch_col = "MS_batch",
                                   par.prior = TRUE,
                                   fill_the_missing = NULL,
                                   no_fit_imputed = TRUE,
                                   qual_col = NULL,
                                   qual_value = NULL,
                                   keep_all = "default",
                                   covariates_cols = NULL) {
    .Deprecated("correct_with_ComBat")
    correct_with_ComBat(
        x = df_long, sample_annotation = sample_annotation,
        feature_id_col = feature_id_col, measure_col = measure_col,
        sample_id_col = sample_id_col, batch_col = batch_col,
        format = "long", par.prior = par.prior,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        keep_all = keep_all, no_fit_imputed = no_fit_imputed,
        qual_col = qual_col, qual_value = qual_value
    )
}

#' @export
#' @rdname correct_batch_effects
correct_with_ComBat_dm <- function(data_matrix, sample_annotation = NULL,
                                   feature_id_col = "peptide_group_label",
                                   measure_col = "Intensity",
                                   sample_id_col = "FullRunName",
                                   batch_col = "MS_batch",
                                   par.prior = TRUE,
                                   fill_the_missing = NULL,
                                   covariates_cols = NULL) {
    .Deprecated("correct_with_ComBat")
    correct_with_ComBat(
        x = data_matrix, sample_annotation = sample_annotation,
        feature_id_col = feature_id_col, measure_col = measure_col,
        sample_id_col = sample_id_col, batch_col = batch_col,
        format = "wide", par.prior = par.prior,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing
    )
}

#' @export
#' @rdname correct_batch_effects
correct_batch_effects_df <- function(df_long, sample_annotation,
                                     continuous_func = NULL,
                                     discrete_func = c(
                                         "MedianCentering", "MeanCentering",
                                         "ComBat", "removeBatchEffect"
                                     ),
                                     batch_col = "MS_batch",
                                     feature_id_col = "peptide_group_label",
                                     sample_id_col = "FullRunName",
                                     measure_col = "Intensity",
                                     order_col = "order",
                                     keep_all = "default",
                                     no_fit_imputed = TRUE,
                                     qual_col = NULL,
                                     qual_value = NULL,
                                     fill_the_missing = NULL,
                                     par.prior = TRUE,
                                     covariates_cols = NULL,
                                     min_measurements = 8, ...) {
    .Deprecated("correct_batch_effects")
    correct_batch_effects(
        x = df_long, sample_annotation = sample_annotation, format = "long",
        continuous_func = continuous_func, discrete_func = discrete_func,
        batch_col = batch_col, feature_id_col = feature_id_col,
        sample_id_col = sample_id_col, measure_col = measure_col,
        order_col = order_col, keep_all = keep_all, no_fit_imputed = no_fit_imputed,
        qual_col = qual_col, qual_value = qual_value,
        fill_the_missing = fill_the_missing, par.prior = par.prior,
        covariates_cols = covariates_cols,
        min_measurements = min_measurements, ...
    )
}

#' @export
#' @rdname correct_batch_effects
correct_batch_effects_dm <- function(data_matrix, sample_annotation,
                                     continuous_func = NULL,
                                     discrete_func = c("MedianCentering", "MeanCentering", "ComBat"),
                                     batch_col = "MS_batch",
                                     feature_id_col = "peptide_group_label",
                                     sample_id_col = "FullRunName",
                                     measure_col = "Intensity",
                                     order_col = "order",
                                     min_measurements = 8,
                                     no_fit_imputed = TRUE,
                                     fill_the_missing = NULL,
                                     par.prior = TRUE,
                                     covariates_cols = NULL,
                                     ...) {
    .Deprecated("correct_batch_effects")
    correct_batch_effects(
        x = data_matrix, sample_annotation = sample_annotation, format = "wide",
        continuous_func = continuous_func, discrete_func = discrete_func,
        batch_col = batch_col, feature_id_col = feature_id_col,
        sample_id_col = sample_id_col, measure_col = measure_col,
        order_col = order_col, min_measurements = min_measurements,
        no_fit_imputed = no_fit_imputed, fill_the_missing = fill_the_missing,
        par.prior = par.prior, covariates_cols = covariates_cols, ...
    )
}

#' @export
#' @rdname correct_batch_effects
correct_with_removeBatchEffect_df <- function(df_long, sample_annotation = NULL,
                                              feature_id_col = "peptide_group_label",
                                              measure_col = "Intensity",
                                              sample_id_col = "FullRunName",
                                              batch_col = "MS_batch",
                                              covariates_cols = NULL,
                                              fill_the_missing = NULL,
                                              keep_all = "default", ...) {
    .Deprecated("correct_with_removeBatchEffect")
    correct_with_removeBatchEffect(
        x = df_long, sample_annotation = sample_annotation,
        feature_id_col = feature_id_col, measure_col = measure_col,
        sample_id_col = sample_id_col, batch_col = batch_col,
        format = "long", covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing, keep_all = keep_all, ...
    )
}

#' @export
#' @rdname correct_batch_effects
correct_with_removeBatchEffect_dm <- function(data_matrix, sample_annotation,
                                              feature_id_col = "peptide_group_label",
                                              measure_col = "Intensity",
                                              sample_id_col = "FullRunName",
                                              batch_col = "MS_batch",
                                              covariates_cols = NULL,
                                              fill_the_missing = NULL, ...) {
    .Deprecated("correct_with_removeBatchEffect")
    correct_with_removeBatchEffect(
        x = data_matrix, sample_annotation = sample_annotation,
        feature_id_col = feature_id_col, measure_col = measure_col,
        sample_id_col = sample_id_col, batch_col = batch_col,
        format = "wide", covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing, ...
    )
}
