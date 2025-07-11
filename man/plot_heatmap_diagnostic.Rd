% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proteome_wide_diagnostics.R
\name{plot_heatmap_diagnostic}
\alias{plot_heatmap_diagnostic}
\title{Plot the heatmap of samples (cols) vs features (rows)}
\usage{
plot_heatmap_diagnostic(
  data_matrix,
  sample_annotation = NULL,
  sample_id_col = "FullRunName",
  factors_to_plot = NULL,
  fill_the_missing = -1,
  color_for_missing = "black",
  heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color_list = NULL,
  peptide_annotation = NULL,
  feature_id_col = "peptide_group_label",
  factors_of_feature_ann = c("KEGG_pathway", "evolutionary_distance"),
  color_list_features = NULL,
  filename = NULL,
  width = 7,
  height = 7,
  units = c("cm", "in", "mm"),
  plot_title = NULL,
  ...
)
}
\arguments{
\item{data_matrix}{features (in rows) vs samples (in columns) matrix, with
feature IDs in rownames and file/sample names as colnames.
See "example_proteome_matrix" for more details (to call the description,
use \code{help("example_proteome_matrix")})}

\item{sample_annotation}{data frame with:
\enumerate{ \item \code{sample_id_col} (this can be repeated as row names)
  \item biological covariates
  \item technical covariates (batches etc) }.
  See \code{help("example_sample_annotation")}}

\item{sample_id_col}{name of the column in \code{sample_annotation} table,
where the filenames (colnames of the \code{data_matrix} are found).}

\item{factors_to_plot}{vector of technical and biological factors to be
plotted in this diagnostic plot (assumed to be present in
\code{sample_annotation})}

\item{fill_the_missing}{numeric value that the missing values are
substituted with, or \code{NULL} if features with missing values are to be
excluded.}

\item{color_for_missing}{special color to make missing values.
Usually black or white, depending on \code{heatmap_color}}

\item{heatmap_color}{vector of colors used in heatmap (typicall a gradient)}

\item{cluster_rows}{boolean value determining if rows
should be clustered}

\item{cluster_cols}{boolean value determining if columns
should be clustered}

\item{color_list}{list, as returned by \code{sample_annotation_to_colors},
where each item contains a color vector for each factor to be mapped to the
color.}

\item{peptide_annotation}{long format data frame with peptide ID and their
corresponding protein and/or gene annotations.
See \code{help("example_peptide_annotation")}.}

\item{feature_id_col}{name of the column with feature/gene/peptide/protein
ID used in the long format representation \code{df_long}. In the wide
formatted representation \code{data_matrix} this corresponds to the row
names.}

\item{factors_of_feature_ann}{vector of factors that characterize features,
as listed in \code{peptide_annotation}}

\item{color_list_features}{list, as returned by
\code{sample_annotation_to_colors},
but mapping \code{peptide_annotation} where each item contains a color vector
for each factor to be mapped to the color.}

\item{filename}{path where the results are saved.
If null the object is returned to the active window;
otherwise, the object is save into the file. Currently only
pdf and png format is supported}

\item{width}{option  determining the output image width}

\item{height}{option  determining the output image height}

\item{units}{units: 'cm', 'in' or 'mm'}

\item{plot_title}{title of the plot (e.g., processing step + representation
level (fragments, transitions, proteins) + purpose (meanplot/corrplot etc))}

\item{...}{other parameters of \code{link[pheatmap]{pheatmap}}}
}
\value{
object returned by \code{link[pheatmap]{pheatmap}}
}
\description{
Plot the heatmap of samples (cols) vs features (rows)
}
\examples{
# Load necessary datasets
data(list = c("example_proteome_matrix", "example_sample_annotation"), package = "proBatch")

log_transformed_matrix <- log_transform_dm(example_proteome_matrix)
heatmap_plot <- plot_heatmap_diagnostic(log_transformed_matrix,
    example_sample_annotation,
    factors_to_plot = c("MS_batch", "digestion_batch", "Diet", "DateTime"),
    cluster_cols = TRUE, cluster_rows = FALSE,
    show_rownames = FALSE, show_colnames = FALSE
)

color_list <- sample_annotation_to_colors(example_sample_annotation,
    factor_columns = c(
        "MS_batch", "EarTag", "Strain",
        "Diet", "digestion_batch", "Sex"
    ),
    numeric_columns = c("DateTime", "order")
)

log_transformed_matrix <- log_transform_dm(example_proteome_matrix)
heatmap_plot <- plot_heatmap_diagnostic(log_transformed_matrix,
    example_sample_annotation,
    factors_to_plot = c("MS_batch", "digestion_batch", "Diet", "DateTime"),
    cluster_cols = TRUE, cluster_rows = FALSE,
    color_list = color_list,
    show_rownames = FALSE, show_colnames = FALSE
)

}
\seealso{
\code{\link{sample_annotation_to_colors}},
\code{\link[pheatmap]{pheatmap}}
}
