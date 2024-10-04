#' Linear Regression R6 Class
#'
#' This R6 class performs linear regression using ordinary least squares.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline stat_summary geom_text labs theme_minimal theme element_blank element_rect element_text element_line unit
#' @importFrom R6 R6Class
#' @field formula A formula object describing the model.
#' @field data A data frame containing the variables in the model.
#' @field X The model matrix generated from the formula and data.
#' @field y The response variable (dependent variable).
#' @field beta_hat The estimated regression coefficients.
#' @field y_hat The fitted values.
#' @field e_hat The residuals.
#' @field df Degrees of freedom of the model.
#' @field sigma_squared Residual variance of the model.
#' @field var_beta_hat Variance of the regression coefficients.
#' @field t_beta t-values for each coefficient.
#'
#' @export
# Define the linreg R6 class
linreg <- R6Class("linreg",
                  public = list(
                    formula = NULL,
                    data = NULL,
                    data_name = NULL,
                    X = NULL,
                    y = NULL,
                    beta_hat = NULL,
                    y_hat = NULL,
                    e_hat = NULL,
                    df = NULL,
                    sigma_squared = NULL,
                    var_beta_hat = NULL,
                    t_beta = NULL,
                    
                    #' @description
                    #' Constructor method for linear regression model
                    #' Performs the linear regression operation
                    #' @param formula A formula object that describes the model
                    #' @param data A data.frame containing the data on which the linear regression is performed
                    #' @return An R6 object representing the linear regression model
                    #' @export
                    initialize = function(formula, data) {
                      self$formula <- formula
                      self$data <- data
                      self$data_name <- deparse(substitute(data))
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
                    
                    #' @description
                    #' Prints the coefficients in an orderly manner
                    print = function() {
                      cat("linreg(formula = ", deparse(self$formula), ", data = ", self$data_name, ")\n", sep = "")
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
                    
                    #' @description
                    #' Returns the values that the model predicted
                                        
                    pred = function() {
                      return(self$y_hat)
                    },
                    
                    #' @description
                    #' Returns the models residuals
                    resid = function() {
                      return(self$e_hat)
                    },
                    
                    #' @description
                    #' Returns the models coefficients as a named vector
                    #' @return coefficients as named vector
                    
                    coef = function() {
                      named_vec <- as.numeric(self$beta_hat)
                      names(named_vec) <- rownames(self$beta_hat)
                      return(named_vec)
                    },
                    
                    #' @description
                    #' Plots the residuals vs the fitted values using ggplot2
                    plot = function() {
                      plot_helper <- function(residuals, plot_title, y_label){
                        # Index the data points to allow for index markers in the plot
                        plot_data <- data.frame(
                          Fitted = self$y_hat,
                          Residuals = residuals,
                          Index = 1:length(self$e_hat)
                        )
                        # Threshold for outliers with large residuals
                        threshold <- 2.5
                        plot_data$Influential <- abs(scale(self$e_hat)) > threshold
                        p <- ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
                          geom_point(shape = 1, color = "black", size = 3, stroke = 1.5) + # Circles 
                          geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 1) + # Dotted y_intercept line
                          stat_summary(fun = median, geom = "line", color = "red", linewidth = 1.2) + # Line intercepting the median
                          geom_text(
                            data = subset(plot_data, Influential == TRUE),
                            aes(label = Index),
                            hjust = +1.3, vjust = -0.3, color = "black", size = 3
                          ) +
                          labs(
                            title = plot_title,
                            x = paste("Fitted values\n lm(", deparse(self$formula),")"),
                            y = y_label
                          ) + 
                          theme_minimal() + 
                          theme(
                            panel.grid.major = element_blank(),   # Remove major grid lines
                            panel.grid.minor = element_blank(),   # Remove minor grid lines
                            panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border around plot
                            axis.text.x = element_text(size = 14),  # Increase x-tick font
                            axis.text.y = element_text(size = 14, angle = 90, vjust = 0.5),  # Increase y-tick font, rotate text
                            axis.title.x = element_text(size = 16), # Increase x-label font size
                            axis.title.y = element_text(size = 16), # Increase y-label font size
                            axis.ticks.x = element_line(linewidth=0.5) ,  # Enable x-ticks
                            axis.ticks.y = element_line(linewidth=0.5),  # Enable y-ticks
                            axis.ticks.length = unit(0.3, "cm"), # Increase tick length
                            plot.title = element_text(size = 20, hjust = 0.5)  # Center the title and increase its size
                          )
                        p <- p + theme(aspect.ratio = 2/3) # Change aspect ratio
                        print(p)
                      }
                      plot_helper(self$e_hat, # The non standardized version
                                  "Residuals", 
                                  "Residuals vs Fitted") 
                      
                      # Standardized sqrt residuals
                      standardized_residuals <- scale(self$resid())
                      sqrt_abs_standardized_residuals <- sqrt(abs(standardized_residuals))
                      
                      plot_helper(
                        residuals = sqrt_abs_standardized_residuals, # Standardize version
                        plot_title = "Scale−Location",
                        y_label = expression(sqrt(abs("Standardized residuals")))
                      )
                    },
                    
                    #' @description
                    #' Placeholder for implementation of summary function
                    summary = function() {
                      cat("TODO: Summary method to be implemented\n")
                    }
                  )
)

