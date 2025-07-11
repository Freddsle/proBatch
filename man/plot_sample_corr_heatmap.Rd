% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation-based_diagnostics.R
\name{plot_sample_corr_heatmap}
\alias{plot_sample_corr_heatmap}
\title{Sample correlation matrix (heatmap)}
\usage{
plot_sample_corr_heatmap(
  data_matrix,
  samples_to_plot = NULL,
  sample_annotation = NULL,
  sample_id_col = "FullRunName",
  factors_to_plot = NULL,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  heatmap_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  color_list = NULL,
  filename = NULL,
  width = NA,
  height = NA,
  units = c("cm", "in", "mm"),
  plot_title = sprintf("Correlation matrix of\%s samples",
    ifelse(is.null(samples_to_plot), "", " selected")),
  ...
)
}
\arguments{
\item{data_matrix}{features (in rows) vs samples (in columns) matrix, with
feature IDs in rownames and file/sample names as colnames.
See "example_proteome_matrix" for more details (to call the description,
use \code{help("example_proteome_matrix")})}

\item{samples_to_plot}{string vector of samples in
\code{data_matrix} to be used in the plot}

\item{sample_annotation}{data frame with:
\enumerate{ \item \code{sample_id_col} (this can be repeated as row names)
  \item biological covariates
  \item technical covariates (batches etc) }.
  See \code{help("example_sample_annotation")}}

\item{sample_id_col}{name of the column in \code{sample_annotation} table,
where the filenames (colnames of the \code{data_matrix} are found).}

\item{factors_to_plot}{vector of technical and biological covariates to be
plotted in this diagnostic plot (assumed to be present in
\code{sample_annotation})}

\item{cluster_rows}{boolean values determining if rows should be clustered or \code{hclust} object}

\item{cluster_cols}{boolean values determining if columns should be clustered or \code{hclust} object}

\item{heatmap_color}{vector of colors used in heatmap.}

\item{color_list}{list, as returned by \code{sample_annotation_to_colors},
where each item contains a color vector for each factor to be mapped to the
color.}

\item{filename}{path where the results are saved.
If null the object is returned to the active window;
otherwise, the object is save into the file. Currently only
pdf and png format is supported}

\item{width}{option  determining the output image width}

\item{height}{option  determining the output image height}

\item{units}{units: 'cm', 'in' or 'mm'}

\item{plot_title}{title of the plot (e.g., processing step + representation
level (fragments, transitions, proteins) + purpose (meanplot/corrplot etc))}

\item{...}{parameters for the \code{\link[pheatmap]{pheatmap}} visualisation, for details see
examples and help to corresponding functions}
}
\value{
\code{pheatmap} object
}
\description{
Plot correlation of selected samples
}
\examples{
data(list = c("example_sample_annotation", "example_proteome_matrix"), package = "proBatch")
specified_samples <- example_sample_annotation$FullRunName[
    which(example_sample_annotation$order \%in\% 110:115)
]

sample_corr_heatmap <- plot_sample_corr_heatmap(example_proteome_matrix,
    samples_to_plot = specified_samples,
    factors_to_plot = c("MS_batch", "Diet", "DateTime", "digestion_batch"),
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_names_col = TRUE, annotation_legend = FALSE,
    show_colnames = FALSE
)

color_list <- sample_annotation_to_colors(example_sample_annotation,
    factor_columns = c(
        "MS_batch", "EarTag", "Strain",
        "Diet", "digestion_batch", "Sex"
    ),
    numeric_columns = c("DateTime", "order")
)
sample_corr_heatmap_annotated <- plot_sample_corr_heatmap(log_transform_dm(example_proteome_matrix),
    sample_annotation = example_sample_annotation,
    factors_to_plot = c("MS_batch", "Diet", "DateTime", "digestion_batch"),
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_names_col = TRUE,
    show_colnames = FALSE, color_list = color_list
)

}
\seealso{
\code{\link[pheatmap]{pheatmap}}
}
