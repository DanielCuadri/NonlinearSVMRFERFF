#' Scale numeric columns of a data frame
#'
#' This function standardizes the numeric columns of a data frame.
#' If means (`m`) and standard deviations (`s`) are provided, it uses them
#' to scale the data. Otherwise, it calculates the mean and standard deviation
#' from the data and returns both the scaled data and the computed statistics.
#'
#' @param data Data frame containing numeric and/or non-numeric variables.
#' @param m (optional) List of pre-computed means for each column.
#' @param s (optional) List of pre-computed standard deviations for each column.
#'
#' @return
#' - If `m` and `s` are provided: returns a scaled data frame.
#' - If not provided: returns a list with:
#'   \item{scaled_data}{Scaled data frame}
#'   \item{mu}{List of means for each column}
#'   \item{sigma}{List of standard deviations for each column}
#'
my_scale <- function(data, m=NULL, s=NULL) {
  
  data <- as.data.frame(data)
  scaled_data <- data
  if (!is.null(m) && !is.null(s)){
    for (col in names(data)) {
      x <- data[[col]]
      
      if (is.numeric(x)) {
        mu <- m[[col]]
        sigma <- s[[col]]
        
        if (sigma == 0) {
          scaled <- x - mu  # only centering
        } else {
          scaled <- (x - mu) / sigma
        }
        
        scaled_data[[col]] <- scaled
      }
    }
    
    return(scaled_data)
  }
  else{
    mu_list <- list()
    sigma_list <- list()
    
    for (col in names(data)) {
      x <- data[[col]]
      
      if (is.numeric(x)) {
        mu <- mean(x, na.rm = TRUE)
        sigma <- sd(x, na.rm = TRUE)
        
        if (sigma == 0) {
          scaled <- x - mu  # only centering
        } else {
          scaled <- (x - mu) / sigma
        }
        
        scaled_data[[col]] <- scaled
        mu_list[[col]] <- mu
        sigma_list[[col]] <- sigma
      }
    }
    
    return(list(
      scaled_data = scaled_data,
      mu = mu_list,
      sigma = sigma_list
    ))
  }
  
}