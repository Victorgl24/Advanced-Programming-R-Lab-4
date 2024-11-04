RidgeRegressor <- R6Class("RidgeRegressor",
                          public = list(
                            formula = NULL,
                            data = NULL,
                            lambda = NULL,
                            coefficients = NULL,
                            fitted_values = NULL,
                            
                            # Initialize method
                            initialize = function(formula, data, lambda) {
                              if (lambda <= 0) stop("Lambda should be a positive number.")
                              self$formula <- formula
                              self$data <- data
                              self$lambda <- lambda
                              self$fit()
                            },
                            
                            # Fit method to calculate ridge regression coefficients
                            fit = function(method = "least_squares") {
                              stopifnot(is.character(method), method %in% c("least_squares", "qr"))
                              if(method == least_squares){
                                X <- model.matrix(self$formula, self$data)
                                y <- model.response(model.frame(self$formula, self$data))
                                
                                # Center and scale X, excluding the intercept (first column)
                                X[, -1] <- scale(X[, -1], center = TRUE, scale = TRUE)
                                
                                # Ridge regression coefficients
                                I <- diag(ncol(X))
                                I[1, 1] <- 0  # Don't regularize the intercept
                                self$coefficients <- solve(t(X) %*% X + self$lambda * I) %*% t(X) %*% y
                                self$fitted_values <- X %*% self$coefficients
                              } else if (method == "qr") {
                                X <- model.matrix(self$formula, self$data)
                                y <- model.response(model.frame(self$formula, self$data))
                                
                                # Center and scale X, excluding the intercept (first column)
                                X[, -1] <- scale(X[, -1], center = TRUE, scale = TRUE)
                                
                                # QR decomposition
                                QR <- qr(X)
                                Q <- qr.Q(QR)
                                R <- qr.R(QR)
                                
                                # Ridge regularization adjustment
                                R_diag <- diag(R)
                                R_diag[-1] <- R_diag[-1] + self$lambda
                                R <- diag(R_diag)
                                
                                # Compute ridge coefficients using QR decomposition
                                self$coefficients <- solve(t(R) %*% R + self$lambda * diag(ncol(R))) %*% t(R) %*% t(Q) %*% y
                                self$fitted_values <- X %*% self$coefficients
                              }
                              
                            },
                            
                            # Print method
                            print = function() {
                              cat("Ridge Regression Coefficients:\n")
                              print(self$coefficients)
                            },
                            
                            # Predict method
                            predict = function(newdata = NULL) {
                              if (is.null(newdata)) {
                                return(self$fitted_values)
                              } else {
                                X_new <- model.matrix(self$formula, newdata)
                                X_new[, -1] <- scale(X_new[, -1], center = TRUE, scale = TRUE)
                                return(X_new %*% self$coefficients)
                              }
                            },
                            
                            # Coef method to return coefficients
                            coef = function() {
                              return(self$coefficients)
                            }
                          )
)
