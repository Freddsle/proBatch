% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation-based_diagnostics.R
\name{calculate_peptide_corr_distr}
\alias{calculate_peptide_corr_distr}
\title{Calculate peptide correlation between and within peptides of one protein}
\usage{
calculate_peptide_corr_distr(
  data_matrix,
  peptide_annotation,
  protein_col = "ProteinName",
  feature_id_col = "peptide_group_label"
)
}
\arguments{
\item{data_matrix}{features (in rows) vs samples (in columns) matrix, with
feature IDs in rownames and file/sample names as colnames.
See "example_proteome_matrix" for more details (to call the description,
use \code{help("example_proteome_matrix")})}

\item{peptide_annotation}{long format data frame with peptide ID and their
corresponding protein and/or gene annotations.
See \code{help("example_peptide_annotation")}.}

\item{protein_col}{column where protein names are specified}

\item{feature_id_col}{name of the column with feature/gene/peptide/protein
ID used in the long format representation \code{df_long}. In the wide
formatted representation \code{data_matrix} this corresponds to the row
names.}
}
\value{
dataframe with peptide correlation coefficients
that are suggested to use for plotting in
\code{\link{plot_peptide_corr_distribution}} as \code{plot_param}:
}
\description{
Calculate peptide correlation between and within peptides of one protein
}
\examples{
data(list = c("example_peptide_annotation", "example_proteome_matrix"), package = "proBatch")
selected_genes <- c("BOVINE_A1ag", "BOVINE_FetuinB", "Cyfip1")
gene_filter <- example_peptide_annotation$Gene \%in\% selected_genes
peptides_ann <- example_peptide_annotation$peptide_group_label
selected_peptides <- peptides_ann[gene_filter]
matrix_test <- example_proteome_matrix[selected_peptides, ]
pep_annotation_sel <- example_peptide_annotation[gene_filter, ]
corr_distribution <- calculate_peptide_corr_distr(matrix_test,
    pep_annotation_sel,
    protein_col = "Gene"
)

}
