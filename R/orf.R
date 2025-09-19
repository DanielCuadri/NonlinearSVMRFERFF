#' Orthogonal Random Fourier Features
#'
#' This function generates orthogonal random Fourier features (ORF) to approximate
#' shift-invariant kernels. ORF improves the variance of standard Random Fourier
#' Features by using orthogonal projection matrices. Supports parallel computation
#' for large feature dimensions.
#'
#' @param X Data frame or matrix of input samples (rows = samples, columns = features).
#' @param D Even integer specifying the total number of random features (must be even).
#' @param gamma Kernel scale parameter (default = 1).
#' @param Xi Optional pre-defined random vectors (matrix). If NULL, orthogonal vectors are generated.
#' @param ncores Number of cores to use for parallel computation when D/2 > n_features.
#'
#' @return A list with:
#'   \item{random_vectors}{Matrix of orthogonal random projection vectors (size D/2 x n_features).}
#'   \item{features}{Transformed feature matrix Z (size n_samples x D).}
#'
#' @details
#' - If D/2 <= n_features, a single orthogonal block is generated via QR decomposition.
#' - If D/2 > n_features, multiple orthogonal blocks are generated in parallel using `ncores`.
#' - Scaling of vectors is applied using chi-squared random variables and `gamma` to match the kernel spectral properties.
#' - Final feature matrix is constructed as \(Z = \sqrt{2/D} [\sin(X W^T), \cos(X W^T)]\).
#'
orthogonal_random_features <- function(X, D, gamma = 1, Xi = NULL, ncores) {
  if (D %% 2 != 0) stop("D must be an even number.")
  
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  half_D <- D / 2
  
  if (is.null(Xi)) {
    if (half_D <= d) {
      
      G <- matrix(rnorm(d * half_D), nrow = d, ncol = half_D)
      Q <- qr.Q(qr(G, LAPACK = TRUE), complete = FALSE)
      W <- t(Q)
    } else {
      # Case D/2 > d: multiple blocks
      block_size <- d
      num_blocks <- ceiling(half_D / block_size)
      
      # Generate all Gaussian blocks at once
      G <- matrix(rnorm(num_blocks * block_size * d), nrow = num_blocks * block_size, ncol = d)
      
      block_idx <- split(seq_len(nrow(G)), rep(1:num_blocks, each = block_size))
      
      cl <- parallel::makeCluster(ncores)
      parallel::clusterExport(cl, varlist = c("G", "block_idx"), envir = environment())
      
      W_list <- parallel::parLapply(cl, block_idx, function(idx) {
        W_block <- G[idx, , drop = FALSE]
        qr.Q(qr(W_block, LAPACK = TRUE), complete = FALSE)
      })
      
      parallel::stopCluster(cl)
      
      W <- do.call(rbind, W_list)
      W <- W[1:half_D, , drop = FALSE]  # truncate to exactly D/2 rows
    }
    
    # SCALING IN PRACTICE
    # W <- sqrt(d) * sqrt(2 * gamma) * W
    
    # SCALING IN THEORY
    r <- sqrt(rchisq(n = nrow(W), df = d))
    W <- sweep(W, 1, r, `*`) * sqrt(2 * gamma)
    colnames(W) <- as.character(1:ncol(W))
    
  } else {
    W <- Xi
  }
  
  # Compute features
  V <- X %*% t(W)
  Z <- sqrt(2 / D) * cbind(sin(V), cos(V))
  
  list(random_vectors = W, features = Z)
}