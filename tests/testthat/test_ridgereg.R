library(MASS)

# Define the test function
test_ridge_regression_fitted_values <- function() {
  # Load data
  data(mtcars)
  
  # Define lambda and formula
  lambda <- 0.1
  formula <- mpg ~ cyl + disp + hp
  
  # Fit ridge regression using MASS::lm.ridge
  ridge_mass <- lm.ridge(formula, data = mtcars, lambda = lambda)
  
  # Calculate fitted values for lm.ridge by multiplying the model matrix with coefficients
  X <- model.matrix(formula, mtcars)
  fitted_values_mass <- X %*% coef(ridge_mass)
  
  # Fit ridge regression using RidgeRegressor (least squares method)
  ridge_ls <- RidgeRegressor$new(formula, mtcars, lambda, method = "least_squares")
  fitted_values_ls <- ridge_ls$predict()  # Get fitted values from RidgeRegressor
  
  # Fit ridge regression using RidgeRegressor (QR decomposition method)
  ridge_qr <- RidgeRegressor$new(formula, mtcars, lambda, method = "qr")
  fitted_values_qr <- ridge_qr$predict()  # Get fitted values from RidgeRegressor
  
  # Set a tolerance for comparison
  tolerance_value <- 1e-2  # Adjust tolerance as needed
  
  # Test that fitted values from least squares method are similar to lm.ridge fitted values
  test_that("Ridge regression fitted values are similar (least squares vs lm.ridge)", {
    expect_equal(as.vector(fitted_values_ls), as.vector(fitted_values_mass), tolerance = tolerance_value)
  })
  
  # Test that fitted values from QR decomposition method are similar to lm.ridge fitted values
  test_that("Ridge regression fitted values are similar (QR vs lm.ridge)", {
    expect_equal(as.vector(fitted_values_qr), as.vector(fitted_values_mass), tolerance = tolerance_value)
  })
  
  # Print message if tests pass
  cat("All fitted values tests passed successfully!\n")
}

# Run the test function
test_ridge_regression_fitted_values()
