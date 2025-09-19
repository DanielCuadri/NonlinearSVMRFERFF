#' Compute Decision Function Values for Linear SVM
#'
#' This function computes the decision function values for a linear SVM
#' in the original feature space, given support vector coefficients and bias.
#'
#' @param X_train Matrix or data frame of training samples (rows = samples, columns = features).
#' @param y_train Vector of training labels (1 or -1).
#' @param alphas Vector of dual coefficients for support vectors.
#' @param b Scalar bias term of the SVM.
#' @param X_val Matrix or data frame of validation/test samples.
#'
#' @return Numeric vector of decision function values for each validation sample.
#'
funcion_decision <- function(X_train, y_train, alphas, b, X_val) {
  # Vector de pesos efectivos en el espacio original
  coef <- alphas * y_train
  
  # Producto escalar de cada val con cada train
  # (X_val %*% t(X_train)) produce matriz n_val × n_train7
  K <- X_val %*% t(X_train)
  
  # Multiplicamos por coef y sumamos b
  f_values <- K %*% coef + b
  
  as.vector(f_values)
}

#' Predict Labels from Decision Function Values
#'
#' This function converts decision function values into class predictions.
#'
#' @param f_values Numeric vector of decision function values.
#'
#' @return Numeric vector of predicted labels (1 or -1).
#'
predecir <- function(f_values) {
  return(ifelse(f_values > 0, 1, -1))
}

#' Evaluate SVM Model Accuracy
#'
#' This function evaluates the classification accuracy of a linear SVM
#' on a validation or test set.
#'
#' @param X_train Matrix or data frame of training samples.
#' @param y_train Vector of training labels (1 or -1).
#' @param alphas Vector of dual coefficients for support vectors.
#' @param b Scalar bias term of the SVM.
#' @param X_val Matrix or data frame of validation/test samples.
#' @param y_val Vector of true labels for validation/test samples.
#'
#' @return Classification accuracy as a percentage (0-100).
#'
evaluar_modelo <- function(X_train, y_train, alphas, b, X_val, y_val) {
  f_values_val <- funcion_decision(X_train, y_train, alphas, b, X_val)
  
  predicciones <- predecir(f_values_val)
  
  precision <- mean(predicciones == y_val) * 100
  return(precision)
}