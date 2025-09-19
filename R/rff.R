#' Random Fourier Features for Kernel Approximation
#'
#' This function generates Random Fourier Features (RFF) to approximate shift-invariant kernels,
#' including Gaussian, Laplace, Cauchy, and Matérn kernels. It returns both the random projection vectors
#' and the transformed feature matrix.
#'
#' @param X Data frame or matrix of input samples (rows = samples, columns = features).
#' @param D Even integer specifying the total number of random features (must be even).
#' @param gamma Kernel scale parameter.
#' @param Xi Optional matrix of pre-defined random vectors. If NULL, random vectors are generated according to the kernel.
#' @param kernel Character string specifying the kernel type: "gaussian", "laplace", "cauchy", or "matern".
#'
#' @return A list with:
#'   \item{random_vectors}{Matrix of random projection vectors (size D/2 x n_features).}
#'   \item{features}{Transformed feature matrix Z (size n_samples x D).}
#'
#' @details
#' - The function uses different distributions to sample random vectors depending on the kernel:
#'   \describe{
#'     \item{Gaussian}{Multivariate normal with variance 2*gamma.}
#'     \item{Laplace}{Cauchy distribution with scale = gamma.}
#'     \item{Cauchy}{Laplace distribution with scale = sqrt(gamma) (implemented manually).}
#'     \item{Matérn}{Multivariate t distribution derived from the Matérn spectral density (nu = 1.5).}
#'   }
#' - The final feature matrix combines sine and cosine projections: \(Z = \sqrt{2/D} [\sin(X \Xi^T), \cos(X \Xi^T)]\).
#' - `D` must be even to split evenly between sine and cosine components.
#'
random_fourier_features <- function(X, D, gamma, Xi = NULL, kernel = "gaussian") {
if (D %% 2 != 0) stop("D must be an even number.")

X_matrix <- as.matrix(X)
n_features <- ncol(X_matrix)

if (is.null(Xi)) {
  if (kernel == "gaussian") {
    xi <- matrix(rnorm((D / 2) * n_features, mean = 0, sd = sqrt(2 * gamma)), nrow = D/2)
    colnames(xi) <- as.character(1:ncol(xi))
  } else if (kernel == "laplace") {
    # Cauchy distribution: location=0, scale=gamma
    xi <- matrix(rcauchy((D / 2) * n_features, location = 0, scale = gamma), nrow = D/2)
    colnames(xi) <- as.character(1:ncol(xi))
  } else if (kernel == "cauchy") {
    # Laplace distribution: location=0, scale=gamma
    # No base R function for Laplace, define it:
    rlaplace <- function(n, location = 0, scale = 1) {
      u <- runif(n, min = -0.5, max = 0.5)
      location - scale * sign(u) * log(1 - 2 * abs(u))
    }
    xi <- matrix(rlaplace((D / 2) * n_features, location = 0, scale = sqrt(gamma)), nrow = D/2)
    colnames(xi) <- as.character(1:ncol(xi))
  } else if (kernel == "matern") {
    nu <- 1.5
    l <- 1 / sqrt(gamma)
    
    n <- nrow(X)
    d <- ncol(X)
    
    alpha <- nu + d/2
    beta <- 2*nu / l^2
    
    Z <- matrix(rnorm(D/2 * d), nrow = D/2, ncol = d)
    U <- rgamma(D/2, shape = alpha, scale = 1/beta)
    xi <- sweep(Z, 1, sqrt(U), FUN = "/")
    
    colnames(xi) <- as.character(1:ncol(xi))
    
  } else {
    stop("Unsupported kernel type.")
  }
} else {
  xi <- Xi
}

projection <- X_matrix %*% t(xi)

sin_part <- sin(projection)
cos_part <- cos(projection)

Z <- sqrt(2 / D) * cbind(sin_part, cos_part)

return(list(random_vectors = xi, features = Z))
}