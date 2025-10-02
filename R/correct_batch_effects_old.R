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
