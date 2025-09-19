#' Recursive Feature Elimination with Random Fourier Features
#'
#' This function performs recursive feature elimination (RFE) using random Fourier
#' feature (RFF) transformations and SVM training with optional kernel types.
#'
#' @param X_train Matrix or data frame of training samples.
#' @param y_train Vector of training labels (numeric: 1/-1 or factor).
#' @param num_runs Number of independent RFE runs (default = 5).
#' @param kernel_choice Character string specifying the kernel type ("linear", "laplace", "ORF" (Gaussian), "cauchy", "IRF"; default = "ORF").
#' @param D_frac Fraction of features s to generate random features: D = s·d (default = 10).
#' @param C_values Numeric vector of candidate SVM cost parameters (default = c(1, 10, 100, 1000, 10000, 1e5)).
#' @param k Number of folds for internal cross-validation (default = 5).
#' @param delF0 Maximum fraction of features eliminated per iteration (default = 0.5).
#' @param delFmin Minimum fraction of features eliminated per iteration (default = 0.01).
#' @param fixed_delF Variable indicating if delF is fixed (1) or iteratively changed from delF0 to delFmin (0) (default = 0).
#' @param seed Random seed for reproducibility (default = 0).
#'
#' @return A list containing:
#'   \item{feature_rankings}{List of feature rankings for each run.}
#'   \item{cv_accuracies}{List of best cross-validation accuracies for each run/iteration.}
#'   
rfe_rff_svm <- function(
    X_train, y_train,
    num_runs = 5,
    kernel_choice = "ORF",
    D_frac = 10,
    C_values = c(1, 10, 100, 1000, 10000, 1e5),
    k = 5,
    delF0 = 0.50,
    delFmin = 0.01,
    fixed_delF = 0,
    seed = 0
) {
  
  # Ensure reproducibility and parallelization
  num_cores <- min(length(C_values), parallel::detectCores())
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  
  resultados <- vector("list", num_runs)
  cv_error <- vector("list", num_runs)
  
  for (i in 1:num_runs) {
    set.seed(seed + i)
    
    # Determine RFF dimension
    if (kernel_choice == "IRF"){
      D <- 4*round((ncol(X_train)*D_frac)/4)
    } else { 
      D <- 2*round((ncol(X_train)*D_frac)/2)
    }
    
    r <- numeric(0)
    mejores_acc <- list()
    s <- 1:ncol(X_train)
    colnames(X_train) <- as.character(1:ncol(X_train))
    transform <- NULL
    sneg <- 0
    spos <- 0
    delF0_p <- delF0
    delFmin_p <- delFmin
    
    # Partition train/validation
    n <- nrow(X_train)
    train_index <- caret::createDataPartition(y_train, p = 0.60, list = FALSE)
    train_y <- as.numeric(as.character(y_train[train_index]))
    val_y   <- as.numeric(as.character(y_train[-train_index]))
    train_x <- X_train[train_index, ]
    val_x   <- X_train[-train_index, ]
    
    while (length(s) > 1) {
      if (length(s) <= 100) delF0_p <- delFmin_p
      delF <- (1 - fixed_delF) * (delFmin_p + (delF0_p - delFmin_p) * (length(s)/ncol(X_train))^0.75) + fixed_delF * delF0_p
      
      data_iter <- train_x[, as.character(s)]
      data_iter2 <- val_x[, as.character(s)]
      sigma2 <- manual_sigest(data_iter, kernel = kernel_choice)
      sigma <- as.numeric((sigma2[1] + sigma2[3]) / 2)
      gamma <- 1 / (2 * sigma^2)
      sigma1 <- sigma
      sigma2_val <- sigma*10
      transform1 <- transform[, as.character(s)]
      
      # RFF / ORF / IRF transformations
      if (kernel_choice == "ORF"){
        result <- orthogonal_random_features(data_iter, D, gamma, transform1, ncores=num_cores)
        transform <- result$random_vectors
        result2 <- orthogonal_random_features(data_iter2, D, gamma, transform1, ncores=num_cores)
      } else if (kernel_choice == "IRF"){
        result <- indefinite_rff(data_iter, D, transform1, sigma1, sigma2_val, sneg=sneg, spos=spos)
        transform <- result$random_vectors
        sneg <- result$sneg
        spos <- result$spos
        result2 <- indefinite_rff(data_iter2, D, transform1, sigma1, sigma2_val, sneg=sneg, spos=spos)
      } else if (kernel_choice != "linear"){
        result <- random_fourier_features(data_iter, D, gamma, transform1, kernel = kernel_choice)
        transform <- result$random_vectors
        result2 <- random_fourier_features(data_iter2, D, gamma, transform1, kernel = kernel_choice)
      } else{
        result <- list(features = as.matrix(data_iter))
        result2 <- list(features = as.matrix(data_iter2))
      }
      
      # Train SVM via svmpath or fallback
      modelo_svm <- tryCatch(
        { svmpath(result$features, train_y) },
        error = function(e) {
          if (grepl("singular", e$message) || grepl("non-conformable", e$message)) return(NULL)
          else stop(e)
        }
      )
      
      if (is.null(modelo_svm)){
        features <- rbind(result$features, result2$features) 
        yy_train <- c(train_y, val_y)
        fold_ids <- sample(rep(1:k, length.out = nrow(features)))
        
        cv_acc <- foreach(C_val = C_values, .combine = c, .packages = "LiblineaR") %dopar% {
          LiblineaR(as.matrix(features), yy_train, type = 1, cost = C_val, cross = k)
        }
        
        best_idx <- which.max(cv_acc)
        best_C <- C_values[best_idx]
        best_acc <- cv_acc[best_idx]*100
        mejores_acc <- append(mejores_acc, best_acc)
        
        capture.output( modelo <- ksvm(x = features, y = yy_train, type = "C-svc", kernel = "vanilladot", C = best_C, scaled = FALSE) ) 
        alphas <- coef(modelo)[[1]] 
        sv_indices <- modelo@SVindex 
        sv <- features[sv_indices, , drop = FALSE] 
        w <- colSums(alphas * sv) 
        w_squared <- w^2
      } else {
        # Evaluate accuracies along SVM path
        accuracies <- numeric(length(modelo_svm$lambda))
        indices_lambda <- round(seq(1, length(modelo_svm$lambda), length.out = min(20, length(modelo_svm$lambda))))
        for (j in indices_lambda){
          coef <- as.vector(modelo_svm$alpha[, j])
          b <- modelo_svm$alpha0[j]
          accuracies[j] <- evaluar_modelo(result$features, train_y, coef, b, result2$features, val_y)
        }
        best_acc <- max(accuracies)
        best_lambda_index <- which.max(accuracies)
        mejores_acc <- append(mejores_acc, best_acc)
        alphas <- as.vector(modelo_svm$alpha[, best_lambda_index])
        v <- alphas * train_y
        w_squared <- as.vector((crossprod(v, result$features))^2)
      }
      
      # Feature elimination scoring function
      if (kernel_choice != "linear"){
        k_idx <- ifelse(1:D > D/2, 1:D - D/2, 1:D)
        abs_xi <- abs(result$random_vectors)
        abs_xi_expanded <- abs_xi[k_idx, , drop = FALSE]
        weighted_matrix <- sweep(abs_xi_expanded, 1, w_squared, `*`)
        t_values <- colSums(weighted_matrix) / D
      } else{
        t_values <- w_squared
      }
      
      dict_t <- setNames(t_values, as.character(s))
      dict_t <- dict_t[order(dict_t)]
      delete <- round(length(s) * delF)
      claves <- names(dict_t[1:delete])
      s <- s[!(s %in% as.numeric(claves))]
      r <- c(r, as.numeric(claves))
    }
    
    r <- c(r, s)
    resultados[[i]] <- r
    cv_error[[i]] <- mejores_acc
  }
  
  parallel::stopCluster(cl)
  
  return(list(feature_rankings = resultados, cv_accuracies = cv_error))
}