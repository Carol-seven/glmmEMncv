#' Penalized Binary GLMM Fit by EM Algorithm
#'
#' @description
#' Fits a binary generalized linear mixed model (GLMM) using an EM algorithm.
#' The fixed-effect coefficients are estimated with penalized regression using
#' lasso, SCAD, or MCP.
#'
#' @param y A binary response vector of length \eqn{n}.
#'
#' @param X An \eqn{n \times p} fixed-effects design matrix.
#'
#' @param W An \eqn{n \times q} random-effects design matrix.
#'
#' @param ni An integer vector specifying the number of observations
#' within each cluster (e.g., subject-specific repeated measurements).
#' The sum of \code{ni} must equal \code{length(y)}.
#'
#' @param Q An integer specifying the number of Gauss–Hermite quadrature points.
#' Larger values improve approximation accuracy but increase computation time.
#'
#' @param penalty A character string (default = "lasso") specifying the penalty.
#' One of \code{"lasso"}, \code{"SCAD"}, or \code{"MCP"}.
#'
#' @param gamma A numeric value (default = NULL) specifying the additional
#' parameter fo the chosen penalty. Defaults to \code{3.7} for SCAD and
#' \code{3} for MCP.
#'
#' @param kfold An integer (default = 5) specifying the number of folds used in
#' cross-validation.
#'
#' @param intercept A logical value (default = FALSE) specifying whether to fit
#' the intercept(s). If \code{TRUE}, an intercept column is added to \code{X}.
#'
#' @param tol A numeric value > 0 (default = 1e-5) specifying the convergence
#' tolerance for the EM algorithm.
#'
#' @param maxiter An integer (default = 1e+3) specifying the maximum number of
#' EM iterations.
#'
#' @param link A character string specifying the link function.
#' Must be one of \code{"logit"} or \code{"probit"}.
#'
#' @return
#' A list containing the following components:
#' \describe{
#' \item{theta}{A vector containing the estimated fixed effects and the
#' random-effect variance.}
#' \item{beta}{Estimated fixed-effects coefficients.}
#' \item{D}{Estimated random-effects covariance matrix.}
#' \item{loglik}{Marginal log-likelihood evaluated by Gauss-Hermite quadrature.}
#' \item{\code{AIC}}{Akaike information criterion.}
#' \item{\code{BIC}}{Bayesian information criterion.}
#' \item{\code{cvfit}}{The final cross-validated \code{ncvreg} model.}
#' \item{iter}{Number of EM iterations.}
#' \item{crit}{Final convergence criterion.}
#' }
#'
#' @details
#' The initial fixed-effect estimates are obtained by applying
#' \code{\link[ncvreg]{cv.ncvreg}} to the binary response using
#' \code{family = "binomial"}. During the EM iterations, the E-step computes
#' the conditional moments of the latent Gaussian variables, and the M-step
#' updates the fixed effects by fitting a penalized Gaussian regression to the
#' working response.
#'
#' For SCAD and MCP, \code{gamma} controls the concavity of the penalty. If not
#' supplied, the default values from common practice are used.
#'
#' @importFrom ncvreg cv.ncvreg
#' @importFrom tmvtnorm mtmvnorm
#'
#' @export

glmmEMncv <- function(y, X, W, ni, Q,
                      penalty = "lasso", gamma = NULL, kfold = 5,
                      intercept = FALSE, tol = 1e-5, maxiter = 1e3,
                      link = c("logit", "probit")) {

  if (is.null(gamma)) {
    gamma.ncv <- switch(penalty, "SCAD" = 3.7, "MCP" = 3, NA)
  } else {
    gamma.ncv <- gamma
  }
  X_full <- if (intercept) cbind(1, X) else X
  p <- ncol(X_full)
  q <- ncol(W)
  m <- length(ni)
  N <- sum(ni)
  cum_sizes <- c(0, cumsum(ni))
  cvfit <- cv.ncvreg(X = X, y = y,
                     family = "binomial", penalty = penalty, gamma = gamma.ncv,
                     nfolds = kfold)
  if (intercept) {
    beta <- matrix(c(as.numeric(cvfit[["fit"]][["beta"]][,cvfit$min])), nrow = p)
  } else {
    beta <- matrix(c(as.numeric(cvfit[["fit"]][["beta"]][-1,cvfit$min])), nrow = p)
  }
  D <- diag(q)
  theta <- c(beta, D[upper.tri(D, diag = TRUE)])
  crit <- Inf
  iter <- 0
  const <- if (link == "logit") 16 * sqrt(3) / (15 * pi) else 1
  while (crit > tol && iter < maxiter) {
    iter <- iter + 1
    sum_b2 <- matrix(0, q, q)
    z_working <- NULL
    for (k in seq_len(m)) {
      idx <- (cum_sizes[k]+1):cum_sizes[k+1]
      y_k <- y[idx]
      X_k <- const * matrix(X_full[idx,], ncol = p)
      W_k <- const * matrix(W[idx,], ncol = q)
      gamma <- X_k %*% beta
      Omega <- W_k %*% D %*% t(W_k) + diag(ni[k])
      Omega_inv <- solve(Omega)
      Delta <- D %*% t(W_k) %*% Omega_inv
      Lambda <- D - D %*% t(W_k) %*% Omega_inv %*% W_k %*% D
      if (length(y_k) == 1) {
        y_k <- matrix(y_k)
      }
      A <- diag(1 - 2 * y_k)
      trunc_sigma <- A %*% Omega %*% A
      trunc_sigma <- (trunc_sigma + t(trunc_sigma)) / 2
      trunc_fit <- mtmvnorm(mean = c(A %*% gamma),
                            sigma = trunc_sigma,
                            upper = rep(0, ni[k]))
      ## E(Z_k | y_k)
      M1 <- A %*% trunc_fit$tmean
      ## E(Z_k Z_k^T | y_k)
      M2 <- A %*% trunc_fit$tvar %*% A + M1 %*% t(M1)
      ## E(b_k | y_k)
      b_k <- Delta %*% (M1 - gamma)
      ## E(b_k b_k^T | y_k)
      b2_k <- Lambda +
        Delta %*%
        (M2 + gamma %*% t(gamma) -
           M1 %*% t(gamma) -
           t(M1 %*% t(gamma))) %*%
        t(Delta)
      sum_b2 <- sum_b2 + b2_k
      z_k <- M1 - W_k %*% b_k
      z_working <- rbind(z_working, z_k)
    }

    theta_old <- theta
    cvfit <- cv.ncvreg(X = X, y = z_working,
                       family = "gaussian", penalty = penalty, gamma = gamma.ncv,
                       nfolds = kfold)
    if (intercept) {
      beta <- matrix(c(as.numeric(cvfit[["fit"]][["beta"]][,cvfit$min])), nrow = p)
    } else {
      beta <- matrix(c(as.numeric(cvfit[["fit"]][["beta"]][-1,cvfit$min])), nrow = p)
    }
    D <- sum_b2 / m
    theta <- c(beta, D[upper.tri(D, diag = TRUE)])
    crit <- sum((theta_old - theta)^2)
  }

  loglik <- glmm_loglik_ghq(y = y, X = X, W = W, beta = beta, ni = ni,
                            Q = Q, link = link)
  npar <- length(theta) - sum(beta == 0)
  AIC <- -2 * loglik + 2 * npar
  BIC <- -2 * loglik + log(N) * npar
  return(list(theta = theta, beta = beta, D = D,
              loglik = loglik, AIC = AIC, BIC = BIC,
              cvfit = cvfit, iter = iter, crit = crit))
}
