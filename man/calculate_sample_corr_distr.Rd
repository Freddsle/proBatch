% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation-based_diagnostics.R
\name{calculate_sample_corr_distr}
\alias{calculate_sample_corr_distr}
\title{Calculates correlation for all pairs of the samples in data matrix, labels
as replicated/same_batch/unrelated in output columns (see "Value").}
\usage{
calculate_sample_corr_distr(
  data_matrix,
  sample_annotation,
  repeated_samples = NULL,
  biospecimen_id_col = "EarTag",
  sample_id_col = "FullRunName",
  batch_col = "MS_batch"
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

\item{repeated_samples}{vector of sample IDs to evaluate, if \code{NULL},
all samples are taken into account for plotting}

\item{biospecimen_id_col}{column in \code{sample_annotation}
that defines a unique bio ID, which is usually a
combination of conditions or groups.
 Tip: if such ID is absent, but can be defined from several columns,
 create new \code{biospecimen_id} column}

\item{sample_id_col}{name of the column in \code{sample_annotation} table,
where the filenames (colnames of the \code{data_matrix} are found).}

\item{batch_col}{column in \code{sample_annotation} that should be used for
batch comparison (or other, non-batch factor to be mapped to color in plots).}
}
\value{
dataframe with the following columns, that
are suggested to use for plotting in
\code{\link{plot_sample_corr_distribution}} as \code{plot_param}:
\enumerate{
\item \code{replicate}
\item \code{batch_the_same}
\item \code{batch_replicate}
\item \code{batches}
}
other columns are: \enumerate{
\item \code{sample_id_1} & \code{sample_id_2}, both
generated from \code{sample_id_col} variable
\item \code{correlation} - correlation of two corresponding samples
\item \code{batch_1} & \code{batch_2} or analogous,
created the same as \code{sample_id_1}
}
}
\description{
Calculates correlation for all pairs of the samples in data matrix, labels
as replicated/same_batch/unrelated in output columns (see "Value").
}
\examples{
data(list = c("example_sample_annotation", "example_proteome_matrix"), package = "proBatch")
corr_distribution <- calculate_sample_corr_distr(
    data_matrix = example_proteome_matrix,
    sample_annotation = example_sample_annotation,
    batch_col = "MS_batch", biospecimen_id_col = "EarTag"
)

}
