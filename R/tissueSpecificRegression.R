#' @title Tissue-specific regression
#'
#' @description Normalizes and performs differential analysis of
#' single-individual transcriptomics data while accouting for
#' relative differences in tissue size between individuals.
#'
#' @param counts RNA-seq biological count \code{matrix} (genes as rows, samples as columns)
#' @param ercc ERCC spike-ins count \code{matrix} (ERCC as rows, samples as columns)
#' @param annots Sample annotations \code{data.frame}
#' @param tissues \code{vector} with tissues as values and genes as names (have to match \code{rownames(counts)})
#' @param mean_design \code{formula} that expresses how \code{counts} depend on variables in \code{annots}
#' @param normalization_groups optional \code{factor} specifying which samples correspond to which group for normalization.
#' Passed to \code{clusters} argument of \link[scran]{computeSumFactors})
#' @param max_iterations \code{integer} specifying maximum number of iterations for maximum-likelihood estimation
#' @param wald_test \code{logical} specifying is Wald test should be performed on regression coefficients
#' @param cores \code{integer} specifying number of cores to be used for parallelization
#' Passed to \code{mc.cores} argument of \link[parallel]{mclapply})
#'
#' @return A list with the following items:
#' \itemize{
#' \item tech_factors: Technical factors (i.e. normalization factors for ERCC-spikeins)
#' \item tissue_size_factors: Tissue-specific size factor (normalization factors for tissue-unique genes divided by normalization factors for ERCC-spikeins)
#' \item final_norm_factors: Gene-specific normalization factors to account for tissue differences.
#' \item coefficients: Regression regressions.
#' \item dispersions: Over-dispersion of negative binomial distribution.
#' \item proportions: Proportion each tissue contributes to total counts of genes.
#' \item convergence: \code{integer} specifying if the model converged? Possible values means that yes. Output of \link[lbfgs]{lbfgs}.
#' \item counts: Raw count \code{matrix}
#' \item norm Normalized count \code{matrix}
#' \item coefficients_se Standard-error of parameter estimates
#' \item pvalue Wald test p-values
#' \item padj Wald test FDR}
#'
#' @details
#' Technical and size factors rely on the normalization procedure implemented as \link[scran]{computeSumFactors}.
#'
#' To estimate regression coefficients, proportions and over-dispersion, we use an iterative approach similar to what is performed in the R function \link[MASS]{glm.nb}. We first get initial estimates of these regression coefficients and proportions assuming that counts are drawn from a Poisson distribution. The iterative procedure consists in estimating the over-dispersion given the mean using the \link[MASS]{theta.ml} function and then in estimating regression coefficients and proportions while keeping the over-dispersion fixed. We repeat this until parameters converge. Optimization is performed using L-BFGS as implemented by \link[lbfgs]{lbfgs}
#'
#' Wald test is performed by using the inverse of the numerical hessian matrix as an estimate of the standard-error.
#'
#' P-values are adjusted using FDR as implemented by \link[stats]{p.adjust}.
#'
#' @example /inst/examples/tissueSpecificRegression.R
#'
#' @export
tissueSpecificRegression <- function(counts,
                                     ercc,
                                     annots,
                                     tissues,
                                     mean_design = ~ 1,
                                     normalization_groups = NULL,
                                     max_iterations = 10L,
                                     wald_test = TRUE,
                                     cores = parallel::detectCores()-1) {

  # compute ercc technical factor
  tech_factors <- normalizationFactor(ercc,
                                      clusters = normalization_groups)

  # compute normalization factors for each tissue
  tissue_norm_factors <- normalizationFactorTissueSpecific(counts,
                                                           tissues = tissues,
                                                           clusters = normalization_groups)
  tissue_norm_factors <- simplify2array(tissue_norm_factors)

  # compute tissue size factor
  tissue_size_factors <- sweep(tissue_norm_factors, 1, tech_factors, "/")

  # create model matrix
  X <- stats::model.matrix(mean_design, annots)
  colnames(X) <- stringr::str_remove_all(colnames(X), "\\(|\\)")
  X <- X[, ! apply(X == 0, 2, all), drop = FALSE]

  # select apply function
  if (cores <= 1) {
    applyfun <- function(X, FUN) pbapply::pblapply(X = X, FUN = FUN)
  } else {
    applyfun <- function(X, FUN) parallel::mclapply(mc.cores = cores, X = X, FUN = FUN)
  }

  # fit model
  fits <-applyfun(seq_len(nrow(counts)), function(i) {

    fitModel(y = counts[i, ],
             X = X,
             tech_factors = tech_factors,
             tissue_size_factors = tissue_size_factors,
             max_iterations = max_iterations,
             wald_test = wald_test)

  })

  # extract results
  convergence <- sapply(fits, function(x) x$convergence)
  names(convergence) <- rownames(counts)

  coefficients <- sapply(fits, function(x) x$beta)
  if (! is.matrix(coefficients)) {
    coefficients <- matrix(coefficients)
    colnames(coefficients) <- colnames(X)
  } else {
    coefficients <- t(coefficients)
  }
  rownames(coefficients) <- rownames(counts)

  proportions <- sapply(fits, function(x) x$prop)
  if (! is.matrix(proportions)) {
    proportions <- matrix(proportions)
    colnames(proportions) <- colnames(tissue_size_factors)
  } else {
    proportions <- t(proportions)
  }
  rownames(proportions) <- rownames(counts)

  delta <- sapply(fits, function(x) x$delta)
  names(delta) <- rownames(counts)

  if (wald_test == TRUE) {

    coefficients_se <- sapply(fits, function(x) x$beta_se)
    if (is.matrix(coefficients_se)) {
      coefficients_se <- t(coefficients_se)
    } else {
      coefficients_se <- matrix(coefficients_se)
    }
    rownames(coefficients_se) <- rownames(counts)
    colnames(coefficients_se) <- colnames(X)

    pvalue <- sapply(fits, function(x) x$pvalue)
    if (is.matrix(pvalue)) {
      pvalue <- t(pvalue)
    } else {
      pvalue <- matrix(pvalue)
    }
    rownames(pvalue) <- rownames(counts)
    colnames(pvalue) <- colnames(X)

    padj <- apply(pvalue, 2, function(x) stats::p.adjust(x, "fdr"))
    rownames(padj) <- rownames(counts)
    colnames(padj) <- colnames(X)

  }

  # final normalization factors
  final_norm_factors <- sweep(proportions %*% t(tissue_size_factors), 2, tech_factors, "*")
  norm <- counts / final_norm_factors

  out <- list(tissue_norm_factors = tissue_norm_factors,
              tech_factors = tech_factors,
              tissue_size_factors = tissue_size_factors,
              final_norm_factors = final_norm_factors,
              coefficients = coefficients,
              dispersions = delta,
              proportions = proportions,
              convergence = convergence,
              counts = counts,
              norm = norm)

  if (wald_test == TRUE) {
    out <- c(out, list(
      coefficients_se = coefficients_se,
      pvalue  = pvalue,
      padj = padj))
  }

  out
}
