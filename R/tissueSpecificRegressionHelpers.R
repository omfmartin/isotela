#' @param params Vector of model parameters: regression coefficients, soft-maxed proportions and inverse overdispersion
#' @param y \code{integer} counts of a gene
#' @param X design \code{matrix} (output of \code{model.matrix(mean_design, annots)})
#' @param log_tech_factors Log of Technical factors (i.e. log of normalization factors for ERCC-spikeins)
#' @param tech_factors Technical factors (i.e. normalization factors for ERCC-spikeins)
#' @param tissue_size_factors Tissue-specific size factor (normalization factors for tissue-unique genes divided by normalization factors for ERCC-spikeins)
#' @param theta Inverse of over-dispersion parameter (size parameter in \link[stats]{dnbinom})
#' @param model Statistical distribution: 'nb' for Negative Binomial or 'poisson' for Poisson
#' @name tissueSpecificRegressionHelpers
NULL

#' @title Model fitting and statistical testing
#'
#' @description Fit negative binomial model using maximum likelihood for a
#' single gene and run Wald test if requested
#'
#' @inheritParams tissueSpecificRegressionHelpers
#' @inheritParams tissueSpecificRegression
#'
#' @return Maximum likelihood estimates of negative binomial model
#' and p-values for regression coefficients if requested
#'
fitModel <- function(y,
                     X,
                     tech_factors,
                     tissue_size_factors,
                     max_iterations = 10,
                     wald_test = TRUE) {

  # initialize parameters
  iter <- 1
  init_converged <- FALSE

  while (init_converged == FALSE & iter <= max_iterations) {
    params_init <- initParams(X, tissue_size_factors)

    init_fit <- lbfgs::lbfgs(vars = params_init,
                             call_eval = computeNegativeLogLikelihood,
                             call_grad = computeNegativeLogLikelihoodGradient,
                             y = y,
                             X = X,
                             tissue_size_factors = tissue_size_factors,
                             log_tech_factors = log(tech_factors),
                             model = "poisson",
                             invisible = 1)
    init_converged <- init_fit$convergence >= 0
    iter <- iter + 1
  }

  # browser()
  # fit negative binomial
  eta <- computeEta(init_fit$par,
                    X = X,
                    log_tech_factors = log(tech_factors),
                    tissue_size_factors = tissue_size_factors)
  theta <- theta_init <- MASS::theta.ml(y = y, mu = exp(eta), limit = 100000)
  # theta <- theta_init <- rgamma(1, 1, 1)
  params_last <- init_fit$par

  epsilon <- 1e-9
  iter <- 1
  change <- Inf

  while(epsilon < change & iter <= max_iterations) {

    fit <- lbfgs::lbfgs(vars = params_last,
                        call_eval = computeNegativeLogLikelihood,
                        call_grad = computeNegativeLogLikelihoodGradient,
                        y = y,
                        X = X,
                        tissue_size_factors = tissue_size_factors,
                        log_tech_factors = log(tech_factors),
                        theta = theta,
                        model = "nb",
                        invisible = 1)

    eta <- computeEta(fit$par,
                      X = X,
                      log_tech_factors = log(tech_factors),
                      tissue_size_factors = tissue_size_factors)
    theta <- MASS::theta.ml(y = y, mu = exp(eta), limit = 10000)
    change <- sum( (params_last - fit$par)^2 )

    params_last <- fit$par
    iter <- iter + 1
  }
  if (iter > max_iterations) {
    fit$convergence <- -99
  }

  # format results
  fit$beta <- fit$par[seq_len(ncol(X))]
  names(fit$beta) <- colnames(X)

  fit$prop <- exp(fit$par[seq.int(ncol(X)+1, length(fit$par))])
  fit$prop <- fit$prop/sum(fit$prop)
  names(fit$prop) <- colnames(tissue_size_factors)

  fit$delta <- as.numeric(1/theta)
  fit$theta <- as.numeric(theta)
  fit$theta_se <- attr(theta, "SE")

  # compute hessian to get standard-error and p-values
  if (wald_test == TRUE) {
    gradient_f <- function(p) {
      computeNegativeLogLikelihood(params = p,
                                   y = y,
                                   log_tech_factors = log(tech_factors),
                                   tissue_size_factors = tissue_size_factors,
                                   model = "nb",
                                   theta = theta,
                                   X = X)
    }

    hessian <- numDeriv::hessian(gradient_f, x = fit$par)

    fit$beta_se <- tryCatch(
      HelpersMG::SEfromHessian(hessian, silent = TRUE)[seq_len(ncol(X))],
      error = function(e) rep_len(NA_real_, ncol(X)))

    fit$pvalue <- stats::pchisq(fit$beta^2/fit$beta_se^2, df = 1, lower.tail = FALSE)
  }

  fit
}

#' @title Linear predictor
#'
#' @description Estimates the mean of the model (i.e. the linear predictor)
#' given the parameter
#'
#' @inheritParams tissueSpecificRegressionHelpers
#'
#' @return Linear predictor
#'
computeEta <- function(params,
                       X,
                       log_tech_factors,
                       tissue_size_factors) {

  index_beta  <- seq_len(ncol(X))
  index_q <- (ncol(X)+1):length(params)

  beta <- params[index_beta]
  prop <- exp(params[index_q])
  prop <- prop/sum(prop)

  eta <- X %*% beta + log_tech_factors + log(tissue_size_factors %*% prop)
  eta
}

#' @title Negative log-likelihood
#'
#' @description Estimates negative log-likelihood of regression
#' coefficients and soft-maxed proportions
#'
#' @inheritParams tissueSpecificRegressionHelpers
#'
#' @title Negative log-likelihoods of model parameters
#'
computeNegativeLogLikelihood <- function(params,
                                         y,
                                         X,
                                         log_tech_factors,
                                         tissue_size_factors,
                                         theta,
                                         model) {


  eta <- computeEta(params,
                    X,
                    log_tech_factors,
                    tissue_size_factors)

  if (model == "nb") {
    ll <- stats::dnbinom(y, mu = exp(eta), size = theta, log = TRUE)
    out <- - sum(ll)

  } else if (model == "poisson") {
    ll <- stats::dpois(y, lambda = exp(eta), log = TRUE)
    out <- - sum(ll)
  }

  out
}

#' @title Negative log-likelihood gradients
#'
#' @description Estimates negative log-likelihood gradients of regression
#' coefficients and soft-maxed proportions
#'
#' @inheritParams tissueSpecificRegressionHelpers
#'
#' @return Negative log-likelihood gradients of model parameters
#'
computeNegativeLogLikelihoodGradient <- function(params,
                                                 y,
                                                 X,
                                                 log_tech_factors,
                                                 tissue_size_factors,
                                                 theta,
                                                 model) {


  index_beta  <- seq_len(ncol(X))
  index_q <- (ncol(X)+1):length(params)

  beta <- params[index_beta]
  exp_q <- exp(params[index_q])
  prop <- exp_q/sum(exp_q)

  eta <- computeEta(params,
                    X,
                    log_tech_factors,
                    tissue_size_factors)
  exp_eta <- exp(eta)

  if (model == "nb") {
    dnll_deta <- (exp_eta - y) / (1 + exp_eta/theta)
  } else if (model == "poisson") {
    dnll_deta <- exp_eta - y
  }

  # gradient beta
  deta_dbeta <- X
  dnll_dbeta <- deta_dbeta * matrix(dnll_deta, nrow(deta_dbeta), ncol(deta_dbeta))
  grad_beta  <- matrixStats::colSums2(dnll_dbeta)

  # gradient q
  sum_phi_exp_q <- tissue_size_factors %*% exp_q
  deta_dq <- tissue_size_factors / matrix(sum_phi_exp_q, length(sum_phi_exp_q), ncol(tissue_size_factors))
  deta_dq <- deta_dq - 1/sum(exp_q)
  deta_dq <- deta_dq * matrix(exp_q, nrow(deta_dq), ncol(deta_dq), byrow = TRUE)
  # colnames(deta_dq) <- colnames(tissue_size_factors)

  # dnll_dq <- sweep(deta_dq, 1, dnll_deta, "*")
  dnll_dq <- deta_dq * matrix(dnll_deta, nrow(deta_dq), ncol(deta_dq))
  grad_q <- matrixStats::colSums2(dnll_dq)

  # return results
  out <- c(grad_beta, grad_q)
  out
}

#' @title  Initialize parameters
#'
#' @description Random starting regression coefficients and
#' soft-maxed proportions used for model fitting
#'
#' @inheritParams tissueSpecificRegressionHelpers
#'
#' @return Random starting parameters
#'
initParams <- function(X, tissue_size_factors) {

  beta_init <- stats::setNames(stats::rnorm(ncol(X)), colnames(X))

  q_init <- log(stats::runif(ncol(tissue_size_factors), min = .3, max = .7))
  q_init <- stats::setNames(q_init, colnames(tissue_size_factors))

  c(beta_init, q_init)
}
