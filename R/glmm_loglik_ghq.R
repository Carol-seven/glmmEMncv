#' GLMM Log-Likelihood via Gauss–Hermite Quadrature
#'
#' @description
#' Computes the marginal log-likelihood for a generalized linear mixed model
#' (GLMM) with binary outcomes using Gauss-Hermite quadrature (GHQ).
#'
#' @param y A binary response vector of length \eqn{n}.
#'
#' @param X An \eqn{n \times p} fixed-effects design matrix.
#'
#' @param W An \eqn{n \times q} random-effects design matrix.
#'
#' @param beta A vector of fixed-effects coefficients of length \eqn{p}.
#'
#' @param ni An integer vector specifying the number of observations
#' within each cluster (e.g., subject-specific repeated measurements).
#' The sum of \code{ni} must equal \code{length(y)}.
#'
#' @param Q An integer specifying the number of Gauss–Hermite quadrature points.
#' Larger values improve approximation accuracy but increase computation time.
#'
#' @param link A character string specifying the link function.
#' Must be one of \code{"logit"} or \code{"probit"}.
#'
#' @return A scalar value representing the marginal log-likelihood of the GLMM.
#'
#' @details
#' The model assumes
#' \deqn{ Y_{ij} \mid b_j \sim \mathrm{Bernoulli}(\pi_{ij}),}
#' where
#' \eqn{g(\pi_{ij}) = X_{ij}^\top \beta + W_{ij}^\top b_j}, and \eqn{g(\cdot)}
#' is either the logit or probit link function.
#'
#' The likelihood contribution for each cluster is obtained by integrating over
#' the random effect distribution using Gauss-Hermite quadrature. The total
#' log-likelihood is the sum of cluster-specific marginal log-likelihoods.
#'
#' @importFrom statmod gauss.quad
#' @importFrom stats plogis pnorm
#'
#' @export

glmm_loglik_ghq <- function(y, X, W, beta, ni,
                            Q, link = c("logit", "probit")) {
  link <- match.arg(link)
  p <- ncol(X)
  q <- ncol(W)
  m <- length(ni)
  ghq <- gauss.quad(Q, kind = "hermite")
  ghq_nodes <- ghq$nodes
  ghq_weights <- ghq$weights
  cluster_loglik <- numeric(m)
  cum_sizes <- c(0, cumsum(ni))
  for (k in seq_len(m)) {
    idx <- (cum_sizes[k]+1):cum_sizes[k+1]
    y_k <- y[idx]
    X_k <- matrix(X[idx,], ncol = p)
    W_k <- matrix(W[idx,], ncol = q)
    eta <- c(X_k %*% beta) + ghq_nodes %*% t(W_k)
    if (link == "probit") {
      log_prob_1 <- pnorm(eta, log.p = TRUE)
      log_prob_0 <- pnorm(-eta, log.p = TRUE)
    } else if (link == "logit") {
      log_prob_1 <- plogis(eta, log.p = TRUE)
      log_prob_0 <- plogis(-eta, log.p = TRUE)
    }
    cond_loglik <- log_prob_1 %*% y_k + log_prob_0 %*% (1 - y_k)
    cluster_loglik[k] <- log(sum(ghq_weights * exp(cond_loglik)))
  }
  return(sum(cluster_loglik))
}
