.post_correction_to_long <- function(corrected_matrix, df_long,
                                     feature_id_col, measure_col, sample_id_col,
                                     original_cols, keep_all) {
    corrected_df <- matrix_to_long(
        corrected_matrix,
        feature_id_col = feature_id_col,
        measure_col    = measure_col,
        sample_id_col  = sample_id_col
    )

    # preserve original values in preBatchCorr_*
    old_measure_col <- .make_pre_col("preBatchCorr", measure_col)
    df_long <- rename(df_long, !!old_measure_col := !!measure_col)

    corrected_df <- left_join(
        corrected_df,
        df_long,
        by = setNames(c(feature_id_col, sample_id_col), c(feature_id_col, sample_id_col))
    )

    default_cols <- c(original_cols, old_measure_col)
    minimal_cols <- c(sample_id_col, feature_id_col, measure_col, old_measure_col)

    subset_keep_cols(
        corrected_df, keep_all,
        default_cols = default_cols,
        minimal_cols = minimal_cols
    )
}
