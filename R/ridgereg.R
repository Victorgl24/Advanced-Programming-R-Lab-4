  RidgeRegressor <- R6Class("RidgeRegressor",
                            public = list(
                              formula = NULL,
                              data = NULL,
                              lambda = NULL,
                              coefficients = NULL,
                              fitted_values = NULL,
                              
                              # Initialize method
                              initialize = function(formula, data, lambda, method = "least_squares"){
                                if (lambda <= 0) stop("Lambda should be a positive number.")
                                self$formula <- formula
                                self$data <- data
                                self$lambda <- lambda
                                self$fit(method)
                              },
                              
                              # Fit method to calculate ridge regression coefficients
                              fit = function(method) {
                                if(method == "least_squares"){
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
  
  library(nycflights13)
  library(dplyr)

  
  AirportDelayVisualizer <- R6Class("AirportDelayVisualizer",
                                    public = list(
                                      airport_data = NULL,
                                      delay_data = NULL,
                                      
                                      # Initialize the class and prepare data
                                      initialize = function() {
                                        self$prepare_data()
                                      },
                                      
                                      # Prepare data: Calculate mean delays and join with airport info
                                      prepare_data = function() {
                                        # Calculate mean delays by destination
                                        self$delay_data <- flights %>%
                                          group_by(dest) %>%
                                          summarise(mean_delay = mean(arr_delay, na.rm = TRUE)) %>%
                                          filter(!is.na(mean_delay))
                                        
                                        # Join with airport locations data
                                        self$airport_data <- self$delay_data %>%
                                          inner_join(airports, by = c("dest" = "faa")) %>%
                                          select(dest, mean_delay, lon, lat)
                                      },
                                      
                                      # Plot the data
                                      plot = function() {
                                        ggplot(self$airport_data, aes(x = lon, y = lat)) +
                                          geom_point(aes(size = mean_delay, color = mean_delay)) +
                                          scale_color_gradient(low = "blue", high = "red") +
                                          labs(
                                            title = "Average Arrival Delay by Airport",
                                            x = "Longitude",
                                            y = "Latitude",
                                            color = "Mean Delay (min)",
                                            size = "Mean Delay (min)"
                                          ) +
                                          theme_minimal()
                                      }
                                    )
  )
  
  visualizer <- AirportDelayVisualizer$new()
  visualizer$plot()
  
   # Load the mtcars dataset
  data(mtcars)
  
  # Initialize the RidgeRegressor with a formula, data, and lambda value
  ridge_model <- RidgeRegressor$new(formula = mpg ~ hp + wt, data = mtcars, lambda = 0.1)
  
  # Print the coefficients
  ridge_model$print()
  
  # Make predictions on the original data
  predicted_values <- ridge_model$predict()
  print(predicted_values)
  
  # Access the coefficients directly
  coefficients <- ridge_model$coef()
  print(coefficients)
