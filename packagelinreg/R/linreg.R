#' Linear Regression R6 Class
#'
#' This R6 class performs linear regression using ordinary least squares.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline stat_summary geom_text labs theme_minimal theme element_blank element_rect element_text element_line unit
#' @importFrom R6 R6Class
#' @field formula A formula object describing the model.
#' @field data A data frame containing the variables in the model.
#' @field data_name contains name of data set
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
linreg <- R6Class("linreg",
                  public = list(
                    # Declare all class attributes as null for code clarity
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
                    #' Constructor method for linear regression model using QR decomposition
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
                      
                      # QR decomposition
                      qr_decomp <- qr(self$X)
                      
                      # Calculate regression coefficients (beta_hat) using qr.coef, avoiding explicit solve()
                      self$beta_hat <- qr.coef(qr_decomp, self$y)
                      
                      # Calculate fitted values (y_hat) and residuals (e_hat)
                      self$y_hat <- self$X %*% self$beta_hat
                      self$e_hat <- self$y - self$y_hat
                      
                      # Degrees of freedom (using rank of the QR decomposition)
                      self$df <- nrow(self$X) - qr_decomp$rank
                      
                      # Residual variance
                      self$sigma_squared <- sum(self$e_hat^2) / self$df
                      
                      # Variance of the regression coefficients
                      R <- qr.R(qr_decomp)  # Upper triangular matrix from QR decomposition
                      self$var_beta_hat <- self$sigma_squared * solve(t(R) %*% R)
                      
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
                    #' Returns the model's residuals
                    resid = function() {
                      return(self$e_hat)
                    },

                    #' @description
                    #' Returns the model's coefficients as a named vector
                    coef = function() {
                      named_vec <- as.numeric(self$beta_hat)
                      names(named_vec) <- rownames(self$beta_hat)
                      return(named_vec)
                    },

                    #' @description
                    #' Plots the residuals vs the fitted values using ggplot2
                    plot = function() {
                      # Helper function is not documented as it is not publically accessible
                      .plot_helper <- function(residuals, plot_title, y_label){
                        plot_data <- data.frame(
                          Fitted = self$y_hat,
                          Residuals = residuals,
                          Index = 1:length(self$e_hat)
                        )
                        threshold <- 2.5 # Threshold for outliers
                        plot_data$Influential <- abs(scale(self$e_hat)) > threshold
                        
                        p <- ggplot(plot_data, aes(x = Fitted, y = Residuals)) +
                          geom_point(shape = 1, color = "black", size = 3, stroke = 1.5) +
                          geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 1) +
                          stat_summary(fun = median, geom = "line", color = "red", linewidth = 1.2) +
                          geom_text(
                            data = subset(plot_data, Influential == TRUE),
                            aes(label = Index),
                            hjust = +1.3, vjust = -0.3, color = "black", size = 3
                          ) +
                          labs( # Set title and axis labels
                            title = plot_title,
                            x = paste("Fitted values\n lm(", deparse(self$formula),")"),
                            y = y_label
                          ) +
                          theme_minimal() +
                          theme(
                            # Remove grid lines, set border, enable ticks, change font sizes
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_rect(color = "black", fill = NA, size = 1),
                            axis.text.x = element_text(size = 14),
                            axis.text.y = element_text(size = 14, angle = 90, vjust = 0.5),
                            axis.title.x = element_text(size = 16),
                            axis.title.y = element_text(size = 16),
                            axis.ticks.x = element_line(linewidth=0.5),
                            axis.ticks.y = element_line(linewidth=0.5),
                            axis.ticks.length = unit(0.3, "cm"),
                            plot.title = element_text(size = 20, hjust = 0.5)
                          )
                        p <- p + theme(aspect.ratio = 2/3) # set aspect ratio
                        print(p)
                      }
                      .plot_helper(self$e_hat, "Residuals", "Residuals vs Fitted") # call helper method using non standardized

                      standardized_residuals <- scale(self$resid())
                      sqrt_abs_standardized_residuals <- sqrt(abs(standardized_residuals))

                      .plot_helper( # standardized
                        residuals = sqrt_abs_standardized_residuals,
                        plot_title = "Scale-Location",
                        y_label = expression(sqrt(abs("Standardized residuals")))
                      )
                    },

                    #' @description
                    #' Summary of the linear regression model
                    summary = function() {
                      coefs <- self$beta_hat
                      std_errors <- sqrt(diag(self$var_beta_hat))
                      t_values <- self$t_beta
                      df <- self$df
                      sigma_squared <- self$sigma_squared
                      residual_se <- sqrt(sigma_squared)

                      p_values <- 2 * (1 - pt(abs(t_values), df = df))

                      significance <- ifelse(p_values < 0.001, "***",
                                             ifelse(p_values < 0.01, "**",
                                                    ifelse(p_values < 0.05, "*", " ")))

                      coef_table <- data.frame(
                        Estimate = format(as.vector(coefs), digits = 4, nsmall = 4),
                        `Std. Error` = format(std_errors, digits = 4, nsmall = 4),
                        `t value` = format(t_values, digits = 4, nsmall = 4),
                        `Pr(>|t|)` = format(p_values, scientific = FALSE, digits = 4, nsmall = 4),
                        Signif = significance,
                        row.names = rownames(coefs)
                      )

                      cat("Coefficients:\n")
                      print(coef_table, row.names = TRUE)

                      cat("\nResidual standard error: ", format(residual_se, digits = 1),
                          " on ", df, " degrees of freedom\n", sep = "")
                    }
                  )
)

data(iris)
model <- linreg$new(Petal.Length ~ Species, iris)
model$plot()
