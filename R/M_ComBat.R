#' The modification to ComBat by Johson, Rabinovic, and Li (Biostatistics, 2007)
#' to shift data to mean and variance of subset of samples rather
#' than grand mean and pooled variance (http://www.biomedcentral.com/1471-2105/16/63).
#' The m-ComBat method was originally developed by Caleb K. Stein in 2014 (https://github.com/SteinCK/M-ComBat)
#' and have been adapted to use inside proBatch package.
#' Reference:
#'   Stein, C.K., Qu, P., Epstein, J. et al. Removing batch effects from
#'   purified plasma cell gene expression microarrays with modified ComBat.
#'   BMC Bioinformatics 16, 63 (2015). https://doi.org/10.1186/s12859-015-0478-3
#'
#'' @title Batch effect correction with modified ComBat (m-ComBat)
#' @description Removes batch-associated linear effects for plotting/unsupervised
#'   tasks. Works for long or wide via \code{format}. Use \code{covariates_cols}
#'   to keep biological effects in the design (not removed).
#' @inheritParams correct_with_ComBat
#' @param covariates_cols optional sample_annotation columns for the design matrix.
#'   These covariates will be preserved during batch correction.
#' @return Matrix if \code{format="wide"}, data.frame if \code{"format="long"} with batch effects removed
#' @seealso \code{\link[limma]{removeBatchEffect}}
#' @examples
#' data(
#'     list = c("example_sample_annotation", "example_proteome_matrix"),
#'     package = "proBatch"
#' )
#' # Wide format
#' corrected_wide <- correct_with_mComBat(
#'     x = example_proteome_matrix,
#'     sample_annotation = example_sample_annotation,
#'     batch_col = "batch",
#'     covariates_cols = c("condition"),
#'     feature_id_col = "peptide_group_label",
#'     format = "wide",
#'     na_action = FALSE
#' )
#' head(corrected_wide)
#' @export
correct_with_mComBat <- function(
    x, sample_annotation = NULL,
    feature_id_col = "peptide_group_label",
    measure_col = "Intensity",
    sample_id_col = "FullRunName",
    batch_col = "MS_batch",
    format = c("long", "wide"),
    covariates_cols = NULL,
    fill_the_missing = NULL,
    keep_all = "default",
    no_fit_imputed = TRUE,
    qual_col = NULL,
    qual_value = NULL,
    mComBat_center = NULL) {
    # stop if mComBat_center is not provided when use_mComBat is TRUE
    if (is.null(mComBat_center)) {
        stop("mComBat_center must be specified when using correct_with_mComBat.")
    }
    # stop if mComBat_center is not present among batches
    if (!is.null(sample_annotation) && !(mComBat_center %in% unique(sample_annotation[[batch_col]]))) {
        stop("mComBat_center must be present in the batch levels of the provided sample_annotation.")
    }

    correct_with_ComBat(
        x,
        sample_annotation = sample_annotation,
        feature_id_col = feature_id_col,
        measure_col = measure_col,
        sample_id_col = sample_id_col,
        batch_col = batch_col,
        format = format,
        covariates_cols = covariates_cols,
        fill_the_missing = fill_the_missing,
        keep_all = keep_all,
        no_fit_imputed = no_fit_imputed,
        qual_col = qual_col,
        qual_value = qual_value,
        mComBat_center = mComBat_center,
        use_mComBat = TRUE
    )
}


################################################################################
# Internal functions
################################################################################
##    M-ComBat
##   Requirements:
##    Load 'sva' package.
##   Input:
##     'dat' = p by n data.frame or matrix , genomic measure matrix
##                 ( dimensions: probe by sample )
##     'batch' = numeric vector of batch association, length n
##     'center' = numeric value of 'gold-standard' batch
##     'mod' = model matrix of potential covariates
## From: https://github.com/SteinCK/M-ComBat
.m_COMBAT <- function(dat, batch, center, mod) {
    batch <- as.factor(batch)
    batchmod <- model.matrix(~ -1 + batch)
    cat("Found", nlevels(batch), "batches\n")
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    }
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    NAs <- any(is.na(dat))
    if (NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"),
            sep = " "
        )
        stop()
    }
    cat("Standardizing Data across genes\n")

    B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))

    # variance of batch of interest
    var.batch <- apply(dat[, batch == center], 1, var)
    var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1 / n.array, n.array)

    grand.mean <- t(n.batches / n.array) %*% B.hat[1:n.batch, ]
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))

    # accounts for covariates here
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }

    # standardized data
    s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, n.array)))

    cat("Fitting L/S model and finding priors\n")
    batch.design <- design[, 1:n.batch]

    gamma.hat <- solve(t(batch.design) %*% batch.design) %*% t(batch.design) %*% t(as.matrix(s.data))

    delta.hat <- NULL
    for (i in batches) {
        delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, na.rm = T))
    }

    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, sva:::aprior)
    b.prior <- apply(delta.hat, 1, sva:::bprior)

    gamma.star <- delta.star <- NULL

    cat("Finding parametric adjustments\n")
    for (i in 1:n.batch) {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], b.prior[i])

        gamma.star <- rbind(gamma.star, temp[1, ])
        delta.star <- rbind(delta.star, temp[2, ])
    }

    cat("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, ] %*% gamma.star)) / (sqrt(delta.star[j, ]) %*% t(rep(1, n.batches[j])))
        j <- j + 1
    }

    bayesdata <- (bayesdata * (sqrt(var.batch) %*% t(rep(1, n.array)))) + matrix(B.hat[center, ], nrow(dat), ncol(dat))

    return(bayesdata)
}
################################################################################
