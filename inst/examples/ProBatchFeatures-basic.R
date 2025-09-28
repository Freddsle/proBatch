# Shared setup for ProBatchFeatures documentation examples -----------------
data("example_ecoli_data", package = "proBatch")

# Extract data
all_metadata <- example_ecoli_data$all_metadata
all_precursors <- example_ecoli_data$all_precursors
all_protein_groups <- example_ecoli_data$all_protein_groups
all_precursor_pg_match <- example_ecoli_data$all_precursor_pg_match

# Keep only essential
rm(example_ecoli_data)

# Construct a ProBatchFeatures object --------------------------------------
pbf <- ProBatchFeatures(
    data_matrix = all_precursors,
    sample_annotation = all_metadata,
    sample_id_col = "FullRunName",
    level = "peptide"
)

# Register a custom step and evaluate it -----------------------------------
pb_register_step("add_one", function(x) x + 1)
pb_eval(pbf, from = "peptide::raw", steps = "add_one")

# Derived objects for downstream helpers -----------------------------------
pbf_logged <- pb_transform(
    pbf,
    from = "peptide::raw",
    steps = c("log2", "medianNorm"),
    store_fast_steps = TRUE
)

# Get information about the object ---------------------------------------
get_operation_log(pbf_logged)
get_chain(pbf_logged)
get_chain(pbf_logged, as_string = TRUE)
pb_pipeline_name(pbf_logged) # the latest pipeline
pb_pipeline_name(pbf_logged, assay = "peptide::raw")

# Access assays and matrices ------------------------------------------------
head(pb_current_assay(pbf_logged)) # the latest assay
head(pb_assay_matrix(pbf_logged)) # the latest matrix
head(pb_assay_matrix(pbf_logged, assay = "peptide::raw")) # the latest matrix
head(pb_as_wide(pbf_logged)) # the latest assay in wide format
head(pb_as_long(pbf_logged)) # the latest assay in long format

# Pipeline evaluation without storing --------------------------------------
pb_eval(
    pbf,
    from = "peptide::raw",
    steps = c("log2", "medianNorm")
)

# Long-format constructor ---------------------------------------------------
long_pbf <- pb_as_long(pbf_logged) # the latest assay in long format

ProBatchFeatures_from_long(
    df_long = long_pbf,
    sample_annotation = all_metadata,
    sample_id_col = "FullRunName",
    feature_id_col = "peptide_group_label",
    level = "peptide"
)


# Aggregate and add levels --------------------------------------------------

pb_aggregate_level(
    pbf,
    from = "peptide::raw",
    feature_var = "ProteinName",
    new_level = "protein"
)

# Add proteins as a new level and link via mapping
#    all_precursor_pg_match has columns: "Precursor.Id", "Protein.Ids"
pbf <- pb_add_level(
    object = pbf,
    from = "peptide::raw",
    new_matrix = all_protein_groups,
    to_level = "protein", # will name "protein::raw" by default
    mapping_df = all_precursor_pg_match,
    from_id = "Precursor.Id",
    to_id = "Protein.Ids",
    map_strategy = "as_is"
)
