#' Linear Regression Function
#'
#' This function performs linear regression using ordinary least squares.
#'
#' @param formula A formula object describing the model.
#' @param data A data frame containing the variables in the model.
#' @return An object of class "linreg" with computed regression statistics.
#'
#' @examples
#' data(mtcars)
#' fit <- linreg(mpg ~ wt + hp, data = mtcars)
#' print(fit)
linreg <- function(formula, data) {

  # Extract model matrix X and dependent variable y
  X <- model.matrix(formula, data)
  y <- data[[all.vars(formula)[1]]]

  # Calculate regression coefficients (beta)
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y

  # Calculate fitted values (y_hat) and residuals (e_hat)
  y_hat <- X %*% beta_hat
  e_hat <- y - y_hat

  # Degrees of freedom
  df <- nrow(X) - ncol(X)

  # Residual variance
  sigma_squared <- sum(e_hat^2) / df

  # Variance of the regression coefficients
  var_beta_hat <- sigma_squared * solve(t(X) %*% X)

  # t-values for each coefficient
  t_beta <- beta_hat / sqrt(diag(var_beta_hat))

  # Create and return the linreg object
  linreg_object <- list(
    coefficients = beta_hat,
    fitted.values = y_hat,
    residuals = e_hat,
    df = df,
    sigma_squared = sigma_squared,
    var_coefficients = var_beta_hat,
    t_values = t_beta
  )

  class(linreg_object) <- "linreg"
  return(linreg_object)
}

#' Print method for linreg objects
#'
#' Prints the coefficients of the model.
#'
#' @param x An object of class "linreg".
#'
print.linreg <- function(x) {
  cat("Coefficients:\n")
  print(x$coefficients)
}
