#' Estimate kernel scale parameter from data
#'
#' This function computes heuristic estimates for the scale parameter
#' (sigma) of a kernel based on the distances between observations.
#' Supports Gaussian and Laplace kernels.
#'
#' @param X Data frame or matrix of observations (rows = samples, columns = features).
#' @param quantiles Numeric vector of quantiles to use for estimating sigma (default: c(0.1, 0.5, 0.9)).
#' @param kernel Character string specifying the kernel type: "gaussian" (L2 norm) or "laplace" (L1 norm).
#'
#' @return Named numeric vector of sigma estimates for the specified quantiles. Names are in the form "Q10", "Q50", etc.
#'
#' @details
#' - For the Laplace kernel, distances are computed using the L1 (Manhattan) norm.
#' - For the Gaussian kernel or other kernel, squared Euclidean distances are used, and the square root is taken after quantile computation.
#' - A minimum threshold is applied to avoid extremely small sigma values (1e-8 for Laplace, 1e-5 for other).
#'
manual_sigest <- function(X, quantiles = c(0.1, 0.5, 0.9), kernel = "gaussian") {
  X <- as.matrix(X)
  if (kernel == "laplace") {
    dists <- as.vector(proxy::dist(X, method = "Manhattan"))
    sigma_estimates <- pmax(quantile(dists, probs = quantiles), 1e-8)
  }
  else {
    dists <- as.vector(dist(X)^2)
    sigma_estimates <- pmax(sqrt(quantile(dists, probs = quantiles)), 1e-5)
  }
  
  names(sigma_estimates) <- paste0("Q", quantiles * 100)
  return(sigma_estimates)
}