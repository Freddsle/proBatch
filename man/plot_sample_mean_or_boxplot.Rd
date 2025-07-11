% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initial_assessment.R
\name{plot_sample_mean_or_boxplot}
\alias{plot_sample_mean_or_boxplot}
\alias{plot_sample_mean}
\alias{plot_boxplot}
\title{Plot per-sample mean or boxplots for initial assessment}
\usage{
plot_sample_mean(
  data_matrix,
  sample_annotation,
  sample_id_col = "FullRunName",
  batch_col = "MS_batch",
  color_by_batch = FALSE,
  color_scheme = "brewer",
  order_col = "order",
  vline_color = "grey",
  facet_col = NULL,
  filename = NULL,
  width = NA,
  height = NA,
  units = c("cm", "in", "mm"),
  plot_title = NULL,
  theme_name = c("classic", "minimal", "bw", "light", "dark"),
  base_size = 20,
  ylimits = NULL
)

plot_boxplot(
  df_long,
  sample_annotation,
  sample_id_col = "FullRunName",
  measure_col = "Intensity",
  batch_col = "MS_batch",
  color_by_batch = TRUE,
  color_scheme = "brewer",
  order_col = "order",
  facet_col = NULL,
  filename = NULL,
  width = NA,
  height = NA,
  units = c("cm", "in", "mm"),
  plot_title = NULL,
  theme_name = c("classic", "minimal", "bw", "light", "dark"),
  base_size = 20,
  ylimits = NULL,
  outliers = TRUE
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

\item{batch_col}{column in \code{sample_annotation} that should be used for
batch comparison (or other, non-batch factor to be mapped to color in plots).}

\item{color_by_batch}{(logical) whether to color points and connecting lines
by batch factor as defined by \code{batch_col}.}

\item{color_scheme}{named vector, names corresponding to unique batch values of
\code{batch_col} in \code{sample_annotation}. Best created with \link{sample_annotation_to_colors}}

\item{order_col}{column in \code{sample_annotation} that determines sample
order. It is used for in initial assessment plots
(\link{plot_sample_mean_or_boxplot}) and  feature-level diagnostics
(\link{feature_level_diagnostics}). Can be `NULL`
   if sample order is irrelevant (e.g. in genomic experiments). For more
   details,
   order definition/inference, see \link{define_sample_order} and
   \link{date_to_sample_order}}

\item{vline_color}{color of vertical lines, typically denoting
different MS batches in ordered runs; should be \code{NULL} for experiments without intrinsic order}

\item{facet_col}{column  in \code{sample_annotation} with a batch factor to
separate plots into facets; usually 2nd to \code{batch_col}. Most meaningful
for multi-instrument MS experiments (where each instrument has its own
order-associated effects (see \code{order_col}) or simultaneous examination
of two batch factors (e.g. preparation day and measurement day).
For single-instrument case should be set to `NULL`}

\item{filename}{path where the results are saved.
If null the object is returned to the active window;
otherwise, the object is save into the file. Currently only
pdf and png format is supported}

\item{width}{option  determining the output image width}

\item{height}{option  determining the output image height}

\item{units}{units: 'cm', 'in' or 'mm'}

\item{plot_title}{title of the plot (e.g., processing step + representation
level (fragments, transitions, proteins) + purpose (meanplot/corrplot etc))}

\item{base_size}{base font size}

\item{ylimits}{range of y-axis to compare two plots side by side, if required.}

\item{df_long}{data frame where each row is a single feature in a single
sample. It minimally has a \code{sample_id_col}, a \code{feature_id_col}
and a \code{measure_col}, but usually also an \code{m_score} (in OpenSWATH
output result file). See \code{help("example_proteome")} for more details.}

\item{measure_col}{if \code{df_long} is among the parameters, it is the
column with expression/abundance/intensity; otherwise, it is used
internally for consistency.}

\item{outliers}{keep (default) or remove the boxplot outliers}
}
\value{
ggplot2 class object. Thus, all aesthetics can be overridden
}
\description{
Plot per-sample mean or boxplots (showing median and quantiles). In ordered samples,
e.g. consecutive MS runs, order-associated effects are visualised.
}
\details{
functions for quick visual assessment of trends associated, overall
  or specific covariate-associated (see \code{batch_col} and \code{facet_col})
}
\examples{
data(list = c(
    "example_proteome", "example_sample_annotation",
    "example_proteome_matrix"
), package = "proBatch")
mean_plot <- plot_sample_mean(example_proteome_matrix, example_sample_annotation,
    order_col = "order", batch_col = "MS_batch"
)

color_list <- sample_annotation_to_colors(example_sample_annotation,
    factor_columns = c("MS_batch"),
    numeric_columns = c("DateTime", "order")
)
plot_sample_mean(example_proteome_matrix, example_sample_annotation,
    order_col = "order", batch_col = "MS_batch", color_by_batch = TRUE,
    color_scheme = color_list[["MS_batch"]]
)

\donttest{
mean_plot <- plot_sample_mean(example_proteome_matrix,
    example_sample_annotation,
    order_col = "order", batch_col = "MS_batch",
    filename = "test_meanplot.png",
    width = 28, height = 18, units = "cm"
)
}

data(list = c("example_proteome", "example_sample_annotation"), package = "proBatch")
boxplot <- plot_boxplot(log_transform_df(example_proteome),
    sample_annotation = example_sample_annotation,
    batch_col = "MS_batch"
)

color_list <- sample_annotation_to_colors(example_sample_annotation,
    factor_columns = c("MS_batch"),
    numeric_columns = c("DateTime", "order")
)
plot_boxplot(log_transform_df(example_proteome),
    sample_annotation = example_sample_annotation,
    batch_col = "MS_batch", color_scheme = color_list[["MS_batch"]]
)

\donttest{
boxplot <- plot_boxplot(log_transform_df(example_proteome),
sample_annotation = example_sample_annotation,
batch_col = "MS_batch", filename = 'test_boxplot.png',
width = 14, height = 9, units = 'in')
}

}
\seealso{
\code{\link[ggplot2]{ggplot}}, \link{date_to_sample_order}
}
