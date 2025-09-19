#' Indefinite RFF
#'
#' This function generates two sets of random vectors, `W_pos` and `W_neg`,
#' for constructing indefinite random Fourier features (RFF). Each set has
#' a different scaling parameter (`sigma1` and `sigma2`) to control their variance.
#'
#' @param D Total number of random features (must be divisible by 4).
#' @param d Dimensionality of input data.
#' @param sigma1 Standard deviation for the positive vectors (default = 1).
#' @param sigma2 Standard deviation for the negative vectors (default = 10).
#'
#' @return A list containing:
#'   \item{W_pos}{Matrix of positive random vectors (d x D/4).}
#'   \item{W_neg}{Matrix of negative random vectors (d x D/4).}
#'   \item{norm_pos}{Normalization factor for positive vectors (default = 1).}
#'   \item{norm_neg}{Normalization factor for negative vectors (default = 1).}
#'
sample_mu_rff <- function(D, d, sigma1 = 1, sigma2 = 10) {
  
  W_pos <- matrix(rnorm((D/4) * d, 0, 1 / sigma1), nrow = d, ncol = D/4)
  W_neg <- matrix(rnorm((D/4) * d, 0, 1 / sigma2), nrow = d, ncol = D/4)
  
  list(W_pos = W_pos, W_neg = W_neg, norm_pos = 1, norm_neg = 1)
}

#' Generate Indefinite Random Fourier Features
#'
#' This function constructs indefinite random Fourier features (RFF) for
#' kernel approximation. It supports combining
#' positive and negative random vectors, complex-like features, and
#' optional pre-defined random vectors.
#'
#' @param X Data frame or matrix of input samples (rows = samples, columns = features).
#' @param D Total number of random features (must be divisible by 4).
#' @param W Optional matrix of pre-defined random vectors. If NULL, vectors are generated internally.
#' @param sigma1 Standard deviation for positive vectors (default = 10000).
#' @param sigma2 Standard deviation for negative vectors (default = 1).
#' @param mode Integer controlling the sampling method (default = 0). If mode != 0, uses `samplingpolysp`.
#' @param sneg Normalization factor for negative vectors (used internally).
#' @param spos Normalization factor for positive vectors (used internally).
#'
#' @return A list containing:
#'   \item{features}{Matrix of transformed features (n_samples x D).}
#'   \item{random_vectors}{Matrix of random vectors used (D x n_features).}
#'   \item{sneg}{Normalization factor for negative vectors.}
#'   \item{spos}{Normalization factor for positive vectors.}
#'
#' @details
#' - If `W` is not provided, random vectors are generated using `sample_mu_rff`
#'   or `samplingpolysp` depending on `mode`.
#' - Combines positive and negative projections to create indefinite kernel features.
#' - Features include cosine and sine transformations scaled according to `spos` and `sneg`.
#'
indefinite_rff <- function(X, D, W = NULL, sigma1 = 10000, sigma2 = 1, mode = 0, sneg = 0, spos = 0) {
  X <- as.matrix(X)
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (is.null(W)) {
    if (mode == 0){
      rff_data <- sample_mu_rff(D, d, sigma1, sigma2)
    }
    else{
      rff_data <- samplingpolysp(D/4, d, p=2, a=2)
    }
    W_pos = rff_data$W_pos
    W_neg = rff_data$W_neg
    spos = rff_data$norm_pos
    sneg = rff_data$norm_neg
    W_pos <- as.matrix(W_pos)
    W_neg <- as.matrix(W_neg)
    
    rownames(W_pos) <- as.character(1:nrow(W_pos))
    rownames(W_neg) <- as.character(1:nrow(W_neg))
  }
  else{
    nc <- nrow(W)
    W_pos <- W[1:(nc/2),]
    W_neg <- W[(nc/2 + 1):nc,]
    W_pos <- t(W_pos)
    W_neg <- t(W_neg)
  }
  
  # Projections
  proj_pos <- X %*% W_pos  # n x D/4
  proj_neg <- X %*% W_neg  # n x D/4
  
  if (sneg < 1e-4){
    Z_pos_real <- sqrt(4 / D) * cos(proj_pos)
    Z_pos_imag <- sqrt(4 * spos / D) * sin(proj_pos)
    
    Z <- cbind(Z_pos_real, Z_pos_imag)
  }
  else{
    Z_pos_real <- sqrt(4 * spos / D) * cos(proj_pos)
    Z_pos_imag <- sqrt(4 * spos / D) * sin(proj_pos)
    
    # For Kernel approximation
    # Z_neg_real <- sqrt(4 * sneg / D) * 1i * cos(proj_neg)
    # Z_neg_imag <- sqrt(4 * sneg / D) * 1i * sin(proj_neg)
    # Z <- cbind(Z_pos_real, Z_pos_imag, -Z_neg_real, -Z_neg_imag)
    
    # For SVM
    Z_neg_real <- sqrt(4 * sneg / D) * cos(proj_neg)
    Z_neg_imag <- sqrt(4 * sneg / D) * sin(proj_neg)
    Z <- cbind(Z_pos_real, Z_pos_imag, Z_neg_real, Z_neg_imag)
  }
  
  W <- cbind(W_pos, W_neg)
  
  return(list(features = Z, random_vectors = t(W), sneg = sneg, spos = spos))
}