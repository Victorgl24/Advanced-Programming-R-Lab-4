#' @title Ridge Regression using last squares and QR decomposition
#' @description
#' This class implements ridge regression using R6 and supports both
#' least-squares and QR decomposition methods for coefficient calculation.
#'
#' @field formula The formula specifying the regression model.
#' @field data The data frame containing the data for the model.
#' @field lambda The regularization parameter; must be a positive number.
#' @field coefficients The computed ridge regression coefficients.
#' @field fitted_values The fitted values from the regression.
#' @importFrom R6 R6Class
#' @export
RidgeRegressor <- R6Class(
  "RidgeRegressor",
  public = list(
    formula = NULL,
    data = NULL,
    lambda = NULL,
    coefficients = NULL,
    fitted_values = NULL,
    
    #' @description Initialize the RidgeRegressor Class
    #' @param formula A formula object specifying the model.
    #' @param data A data frame containing the variables in the formula.
    #' @param lambda A positive numeric value for regularization strength.
    #' @param method The method for calculating coefficients, either
    #'        least_squares (default) or qr for QR decomposition.
    #' @return An object of class RidgeRegressor.
    #' @examples
    #' data(mtcars)
    #' model <- RidgeRegressor$new(formula = mpg ~ hp + wt, data = mtcars, lambda = 0.1)
    #' @export
    initialize = function(formula, data, lambda, method = "least_squares") {
      if (lambda <= 0)
        stop("Lambda should be a positive number.")
      self$formula <- formula
      self$data <- data
      self$lambda <- lambda
      self$fit(method)
    },
    
    #' @description Fit Method for Ridge Regression
    #' Calculates ridge regression coefficients using either least-squares
    #' or QR decomposition based on the specified method.
    #' @param method Character string indicating the fitting method to use. Choices are \"least_squares\" or \"qr\".
    fit = function(method) {
      if (method == "least_squares") {
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
    
    #' Print Method for Coefficients
    #'
    #' Displays the ridge regression coefficients.
    print = function() {
      cat("Ridge Regression Coefficients:\n")
      print(self$coefficients)
    },
    
    #' Predict Method
    #'
    #' Predicts the response variable for new data using the fitted model.
    #'
    #' @param newdata A data frame with the same predictors as the training data.
    #' @return A numeric vector of predictions.
    #' @examples
    #' new_data <- data.frame(hp = c(110, 120), wt = c(2.5, 3.0))
    #' predictions <- model$predict(newdata = new_data)
    predict = function(newdata = NULL) {
      if (is.null(newdata)) {
        return(self$fitted_values)
      } else {
        X_new <- model.matrix(self$formula, newdata)
        X_new[, -1] <- scale(X_new[, -1], center = TRUE, scale = TRUE)
        return(X_new %*% self$coefficients)
      }
    },
    
    #' @description Coef Method
    #' Returns the ridge regression coefficients.
    #' @return A numeric vector of coefficients.
    coef = function() {
      return(self$coefficients)
    }
  )
)

  
#' @importFrom dplyr group_by summarise filter inner_join select %>%
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient labs theme_minimal
#' @export
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
                                      self$delay_data <- nycflights13::flights %>%
                                        group_by(dest) %>%
                                        summarise(mean_delay = mean(arr_delay, na.rm = TRUE)) %>%
                                        filter(!is.na(mean_delay))
                                      
                                      # Join with airport locations data
                                      self$airport_data <- self$delay_data %>%
                                        inner_join(nycflights13::airports, by = c("dest" = "faa")) %>%
                                        select(dest, mean_delay, lon, lat)
                                    },
                                    
                                    # Plot the data
                                    plot = function() {
                                      p <- ggplot(self$airport_data, aes(x = lon, y = lat)) +
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
                                      
                                      print(p)
                                    }
                                  )
)

visualizer <- AirportDelayVisualizer$new()
visualizer$plot()
