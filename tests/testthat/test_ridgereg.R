
# Define the test function
test_ridge_regression <- function() {
  # Load data
  data(mtcars)
  
  # Define lambda and formula
  lambda <- 0.1
  formula <- mpg ~ cyl + disp + hp
  
  # Fit ridge regression using MASS::lm.ridge
  ridge_mass <- lm.ridge(formula, data = mtcars, lambda = lambda)
  
  # Fit ridge regression using RidgeRegressor (least squares method)
  ridge_ls <- RidgeRegressor$new(formula, mtcars, lambda, method = "least_squares")
  
  # Fit ridge regression using RidgeRegressor (QR decomposition method)
  ridge_qr <- RidgeRegressor$new(formula, mtcars, lambda, method = "qr")
  
  # Test that coefficients are similar
  test_that("Ridge regression coefficients are similar (least squares vs lm.ridge)", {
    expect_equal(as.vector(ridge_ls$coef()), as.vector(coef(ridge_mass)), tolerance = 1e-4)
  })
  
  test_that("Ridge regression coefficients are similar (QR vs lm.ridge)", {
    expect_equal(as.vector(ridge_qr$coef()), as.vector(coef(ridge_mass)), tolerance = 1e-4)
  })
  
  # Print messages if tests pass
  cat("All tests passed successfully!\n")
}

# Run the test function
test_ridge_regression()
