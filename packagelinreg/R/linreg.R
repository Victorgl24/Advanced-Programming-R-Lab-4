# Load the R6 package
library(R6)
library(ggplot2)

# Define the linreg R6 class
LinReg <- R6Class("LinReg",
                  public = list(
                    formula = NULL,
                    data = NULL,
                    X = NULL,
                    y = NULL,
                    beta_hat = NULL,
                    y_hat = NULL,
                    e_hat = NULL,
                    df = NULL,
                    sigma_squared = NULL,
                    var_beta_hat = NULL,
                    t_beta = NULL,
                    
                    # Perform linear regression
                    initialize = function(formula, data) {
                      self$formula <- formula
                      self$data <- data
                      self$X <- model.matrix(formula, data)
                      self$y <- data[[all.vars(formula)[1]]]
                      
                      # Calculate regression coefficients (beta)
                      self$beta_hat <- solve(t(self$X) %*% self$X) %*% t(self$X) %*% self$y
                      
                      # Calculate fitted values (y_hat) and residuals (e_hat)
                      self$y_hat <- self$X %*% self$beta_hat
                      self$e_hat <- self$y - self$y_hat
                      
                      # Degrees of freedom
                      self$df <- nrow(self$X) - ncol(self$X)
                      
                      # Residual variance
                      self$sigma_squared <- sum(self$e_hat^2) / self$df
                      
                      # Variance of the regression coefficients
                      self$var_beta_hat <- self$sigma_squared * solve(t(self$X) %*% self$X)
                      
                      # t-values for each coefficient
                      self$t_beta <- self$beta_hat / sqrt(diag(self$var_beta_hat))
                    },
                    
                    # Print method for linreg object
                    print = function() {
                      coef_names <- rownames(self$beta_hat)
                      values <- as.vector(self$beta_hat)
                      max_width <- max(nchar(coef_names), nchar(format(values, digits = 3, nsmall = 2)))
                      
                      cat("Coefficients:\n")
                      for (name in coef_names) {
                        cat(format(name, width = max_width, justify = "right"), " ")
                      }
                      cat("\n")
                      
                      for (value in values) {
                        cat(format(value, digits = 3, nsmall = 2, width = max_width, justify = "right"), " ")
                      }
                      cat("\n")
                    },
                    
                    # Getter method for predicted values
                    predict = function() {
                      return(self$y_hat)
                    },
                    
                    # Getter method for residuals
                    residuals = function() {
                      return(self$e_hat)
                    },
                    
                    # Getter method for coefficients
                    coefficients = function() {
                      named_vec <- as.numeric(self$beta_hat)
                      names(named_vec) <- rownames(self$beta_hat)
                      return(named_vec)
                    },
                    
                    # Plot method for linreg object
                    plot = function() {
                      plot_data <- data.frame(
                        Fitted = self$y_hat,
                        Residuals = self$e_hat,
                        Index = 1:length(self$e_hat)
                      )
                      
                      # Threshold for outliers with large residuals
                      threshold <- 2.5
                      plot_data$Influential <- abs(scale(self$e_hat)) > threshold
                      
                      p <- ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
                        geom_point(shape = 1, color = "black", size = 3, stroke = 1.5) +  # Open circle points
                        geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 1) +
                        stat_summary(fun = median, geom = "line", aes(group = 1), color = "red", linewidth = 1.2) +
                        geom_text(
                          data = subset(plot_data, Influential == TRUE),
                          aes(label = Index),
                          hjust = +1.3, vjust = -0.3, color = "black", size = 3
                        ) +
                        labs(
                          title = "Residuals vs Fitted",
                          x = paste("Fitted values\n lm(", deparse(self$formula),")"),
                          y = "Residuals"
                        ) + 
                        theme_minimal() + 
                        theme(
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_rect(color = "black", fill = NA, size = 1),
                          axis.text.y = element_text(angle = 90, vjust = 0.5)
                        )
                      
                      print(p)
                    },
                    
                    # Summary method placeholder
                    summary = function() {
                      cat("TODO: Summary method to be implemented\n")
                    }
                  )
)

# Example usage
data(iris)
mod_object <- LinReg$new(Petal.Length ~ Species, data = iris)
mod_object$print()
mod_object$plot()
