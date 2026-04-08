# ============================================================
# Helper Function: Bootstrap BDS Test Function
# ============================================================
bds_test_bootstrap <- function(residuals, m = 2, eps = 1, n_bootstrap = 499) {
  
  # Standardize residuals
  residuals_std <- (residuals - mean(residuals)) / sd(residuals)
  
  # Compute BDS statistic on actual residuals
  bds_actual <- bds.test(residuals_std, m = m, eps = eps)
  bds_actual_stat <- as.numeric(bds_actual$statistic)
  
  # Initialize bootstrap statistics
  bds_bootstrap_stats <- numeric(n_bootstrap)
  
  # Generate bootstrap resamples
  for (i in 1:n_bootstrap) {
    # Resample with replacement
    idx_bootstrap <- sample(1:length(residuals_std), 
                            size = length(residuals_std), 
                            replace = TRUE)
    residuals_bootstrap <- residuals_std[idx_bootstrap]
    
    # Compute BDS on bootstrap resample
    tryCatch({
      bds_boot <- bds.test(residuals_bootstrap, m = m, eps = eps)
      bds_bootstrap_stats[i] <- as.numeric(bds_boot$statistic)
    }, error = function(e) {
      bds_bootstrap_stats[i] <<- NA
    })
  }
  
  # Calculate empirical p-value
  # Two-tailed: proportion of bootstrap statistics >= actual statistic
  p_value_empirical <- mean(abs(bds_bootstrap_stats) >= abs(bds_actual_stat), na.rm = TRUE)
  
  # Return results
  return(list(
    bds_statistic = bds_actual_stat,
    bds_pvalue_parametric = as.numeric(bds_actual$p.value),  # Original p-value
    bds_pvalue_empirical = p_value_empirical,                # Bootstrap p-value
    bootstrap_stats = bds_bootstrap_stats,
    n_bootstrap = n_bootstrap,
    m = m,
    eps = eps
  ))
}
