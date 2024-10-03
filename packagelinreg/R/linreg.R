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
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y # beta_hat = (X^T*X)-1*X^t*y 

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
print.linreg <- function(object) {
  # Extract coefficient names and values
  coefs <- object$coefficients
  coef_names <- rownames(coefs)
  values <- as.vector(coefs)
  
  max_width <- max(nchar(coef_names), nchar(format(values, digits = 3, nsmall = 2)))
  
  # Format the output
  cat("Coefficients:\n")
  
  for (name in coef_names) {
    cat(format(name, width = max_width, justify = "right"), " ")
  }
  cat("\n")
  
  # Print the values in a single row, right-aligned under the corresponding names
  for (value in values) {
    cat(format(value, digits = 3, nsmall = 2, width = max_width, justify = "right"), " ")
  }
  cat("\n")
}

pred.linreg <- function(object){
  return(object$fitted_values)
}

resid.linreg <- function(object){
  return (object$residuals)
}

coef.linreg <- function(object){
  named_vec = c(as.numeric(object$coefficients))
  names(named_vec) = rownames(object$coefficients)
  return (named_vec)
}

#' Plot method for linreg objects using ggplot2
#'
#' Plots the residuals vs fitted values using ggplot2, and labels influential points.
#'
#' @param object An object of class "linreg".
#' @return A ggplot2 residuals vs fitted plot with annotated influential points.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth stat_summary geom_text labs theme_minimal theme element_text
#' @export
plot.linreg <- function(object) {
  # Create a data frame for ggplot
  plot_data <- data.frame(
    Fitted = object$fitted.values,
    Residuals = object$residuals,
    Index = 1:length(object$residuals)  # Add the observation index
  )
  
  # Set a threshold for labeling large residuals (e.g., standardized residuals > 2)
  threshold <- 2.685
  plot_data$Influential <- abs(scale(object$residuals)) > threshold
  
  # Create the Residuals vs Fitted plot using ggplot2
  p <- ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
    geom_point(shape = 1, color = "black", size = 3, stroke = 1.5) +  # Open circle points
    
    # Add a horizontal line at 0
    geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 1) +
    
    # Add segments to connect medians of residuals for each fitted value
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "red", size = 1.2) +
    
    # Annotate influential points with their observation index
    ggplot2::geom_text(
      data = subset(plot_data, Influential == TRUE),  # Only label influential points
      aes(label = Index),
      hjust = +1.3, vjust = -0.3, color = "black", size = 3
    ) +
    
    labs(
      title = "Residuals vs Fitted",
      x = "Fitted values\n lm(Petal.Length~Species)",
      y = "Residuals"
    ) + 
    theme_minimal() + 
    theme(
      panel.grid.major = element_blank(),   # Remove major grid lines
      panel.grid.minor = element_blank(),   # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border around plot
      axis.text.y = element_text(angle = 90, vjust = 0.5)  # Rotate y-axis labels 90 degrees
    )
  print(p)
}



summary.linreg <- function(object){
  # TODO: implement function
}
data(iris)

mod_object <- linreg(Petal.Length~Species, data = iris)
plot.linreg(mod_object)



