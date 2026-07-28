.post_correction_to_long <- function(corrected_matrix, df_long,
                                     feature_id_col, measure_col, sample_id_col,
                                     original_cols, keep_all,
                                     prefix = "preBatchCorr",
                                     join_method = c("left_join", "merge")) {
    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    # preserve original values in pre*
    old_measure_col <- .make_pre_col(prefix, measure_col)
    df_long <- rename(df_long, !!old_measure_col := !!measure_col)

    join_method <- match.arg(join_method)
    if (identical(join_method, "merge")) {
        corrected_df <- merge(
            corrected_df,
            df_long,
            by = c(feature_id_col, sample_id_col)
        )
    } else {
        corrected_df <- left_join(
            corrected_df,
            df_long,
            by = setNames(c(feature_id_col, sample_id_col), c(feature_id_col, sample_id_col))
        )
    }

    default_cols <- c(original_cols, old_measure_col)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    subset_keep_cols(
        corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}
