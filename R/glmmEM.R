#' Binary GLMM Fit by EM Algorithm
#'
#' @description
#' Fits a binary generalized linear mixed model (GLMM) using an EM algorithm.
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
#' \item{se}{Approximate standard errors for \code{theta}.}
#' \item{loglik}{Marginal log-likelihood evaluated by Gauss-Hermite quadrature.}
#' \item{iter}{Number of EM iterations.}
#' \item{crit}{Final convergence criterion.}
#' }
#'
#' @details
#' The EM algorithm uses a latent-normal representation of the binary GLMM.
#' In the E-step, conditional moments of the latent response and random effects
#' are computed using moments of a truncated multivariate normal distribution.
#' In the M-step, the fixed-effects coefficients and random-effect variance are
#' updated from these conditional moments.
#'
#' This implementation currently supports only one random effect, such as a
#' random intercept model.
#'
#' @importFrom tmvtnorm mtmvnorm
#'
#' @export

glmmEM <- function(y, X, W, ni, Q,
                   tol = 1e-5, maxiter = 1e3,
                   link = c("logit", "probit")) {
  p <- ncol(X)
  q <- ncol(W)
  m <- length(ni)
  cum_sizes <- c(0, cumsum(ni))
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  D <- diag(q)
  theta <- c(beta, D[upper.tri(D, diag = TRUE)])
  const <- if (link == "logit")  16 * sqrt(3) / (15 * pi) else 1
  crit <- Inf
  iter <- 0
  while (crit > tol && iter < maxiter) {
    iter <- iter + 1
    sum_XtX <- matrix(0, p, p)
    sum_Xtz <- matrix(0, p, 1)
    sum_b2 <- matrix(0, q, q)
    info_mat <- matrix(0, p+1, p+1)
    score <- matrix(0, p+1, 1)

    for (k in seq_len(m)) {
      idx <- (cum_sizes[k]+1):cum_sizes[k+1]
      y_k <- y[idx]
      X_k <- const * matrix(X[idx,], ncol = p)
      W_k <- const * matrix(W[idx,], ncol = q)
      gamma <- X_k %*% beta
      Omega <- W_k %*% D %*% t(W_k) + diag(ni[k])
      Omega_inv <- solve(Omega)
      D_inv <- solve(D)
      Delta <- D %*% t(W_k) %*% Omega_inv
      Lambda <- D - D %*% t(W_k) %*% Omega_inv %*% W_k %*% D
      if (length(y_k) == 1) {
        y_k <- matrix(y_k)
      }
      A <- diag(1 - 2*y_k)
      trunc_sigma <- A %*% Omega %*% A
      trunc_sigma <- (trunc_sigma + t(trunc_sigma)) / 2
      trunc_fit <- mtmvnorm(mean = c(A %*% gamma),
                            sigma = trunc_sigma,
                            upper = rep(0, ni[k]))
      ## E(Z_k \mid y_k)
      M1 <- A %*% trunc_fit$tmean
      ## E(Z_k Z_k^\top \mid y_k)
      M2 <- A %*% trunc_fit$tvar %*% A + M1 %*% t(M1)
      ## E(b_k \mid y_k)
      b_k <- Delta %*% (M1 - gamma)
      ## E(b_k b_k^\top \mid y_k)
      b2_k <- Lambda + Delta %*% (M2 + gamma %*% t(gamma) - M1 %*% t(gamma) - t(M1 %*% t(gamma))) %*% t(Delta)
      sum_XtX <- sum_XtX + t(X_k) %*% X_k
      sum_Xtz <- sum_Xtz + t(X_k) %*% (M1 - W_k %*% b_k)
      sum_b2 <- sum_b2 + b2_k
      score[1:p] <- t(M1 - W_k %*% b_k - gamma) %*% X_k
      score[p+1] <- -0.5 * sum(diag(D_inv %*% matrix(1, 1, 1))) +
        0.5 * sum(diag(D_inv %*% matrix(1, 1, 1) %*% D_inv %*% b2_k))
      info_mat <- info_mat + score %*% t(score)
    }
    theta_old <- theta
    beta <- solve(sum_XtX) %*% sum_Xtz
    D <- sum_b2 / m
    theta <- c(beta, D[upper.tri(D, diag = TRUE)])
    crit <- sum((theta_old - theta)^2)
  }
  loglik <- glmm_loglik_ghq(y, X, W, beta, ni, Q, link)
  se <- sqrt(diag(solve(info_mat)))
  return(list(theta = theta, beta = beta, D = D,
              se = se, loglik = loglik, iter = iter, crit = crit))
}
