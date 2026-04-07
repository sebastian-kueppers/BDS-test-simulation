############################################################
## VAR DATA SIMULATION UTILITIES
## Author: ---
## Purpose: Simulate multivariate VAR(1) data and
##          discretize to Likert-type scales
############################################################

## =========================================================
## 1. Min–max scaling to Likert scale
## =========================================================

min_max_to_scale <- function(x, new_min = 1, new_max = 7) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  # Edge case: constant time series
  if (x_max == x_min) {
    mid <- round((new_min + new_max) / 2)
    return(as.integer(rep(mid, length(x))))
  }
  
  # Linear rescaling
  x_scaled <- new_min + (x - x_min) *
    ((new_max - new_min) / (x_max - x_min))
  
  # Round and clip
  x_round <- pmin(new_max, pmax(new_min, round(x_scaled)))
  
  as.integer(x_round)
}


## =========================================================
## 2. Simulate a single time series
## =========================================================

simulate_ar1 <- function(
    T,
    phi,        # AR(1) coefficient (scalar)
    sigma = 1,  # innovation SD
    burnin = 0,
    init = NULL
) {
  stopifnot(
    is.numeric(phi), length(phi) == 1,
    abs(phi) < 1,   # stationarity condition
    sigma > 0
  )
  
  TT <- T + burnin
  Y  <- numeric(TT)
  
  Y[1] <- if (is.null(init)) rnorm(1) else init
  
  for (t in 2:TT) {
    Y[t] <- phi * Y[t - 1] + rnorm(1, mean = 0, sd = sigma)
  }
  
  if (burnin > 0) Y <- Y[(burnin + 1):TT]
  
  data.frame(V1 = Y)
}


## =========================================================
## 3. Simulate data for N subjects
## =========================================================

simulate_ar1_subjects <- function(
    N,
    T,
    phi,
    sigma = 1,
    burnin = 0
) {
  lapply(seq_len(N), function(i) {
    simulate_ar1(
      T      = T,
      phi    = phi,
      sigma  = sigma,
      burnin = burnin
    )
  })
}



## =========================================================
## 4. Discretize subject data to Likert scales
## =========================================================

discretize_likert <- function(
    raw_data,
    scale_min = 1,
    scale_max = 7,
    suffix = NULL
) {
  likert_data <- lapply(raw_data, function(df) {
    out <- as.data.frame(
      lapply(df, min_max_to_scale,
             new_min = scale_min,
             new_max = scale_max)
    )
    if (!is.null(suffix)) {
      colnames(out) <- paste0(colnames(out), suffix)
    }
    out
  })
  
  likert_data
}


## =========================================================
## 5. Full simulation grid: N × T
## =========================================================

simulate_ar1_grid <- function(
    N_list,
    T_list,
    phi,
    sigma = 1,
    burnin = 0,
    likert_scales = list(
      "likert_1_7"   = c(1, 7),
      "likert_1_100" = c(1, 100)
    ),
    seed    = NULL,
    verbose = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  data_list <- list()
  
  for (N in N_list) {
    data_list[[paste0("N", N)]] <- list()
    
    for (T in T_list) {
      
      raw_data <- simulate_ar1_subjects(
        N      = N,
        T      = T,
        phi    = phi,
        sigma  = sigma,
        burnin = burnin
      )
      
      # Discretise to Likert scales (reuses your existing function)
      likert_data <- list()
      for (nm in names(likert_scales)) {
        sc <- likert_scales[[nm]]
        likert_data[[nm]] <- discretize_likert(
          raw_data,
          scale_min = sc[1],
          scale_max = sc[2],
          suffix    = paste0("_", nm)
        )
      }
      
      data_list[[paste0("N", N)]][[paste0("T", T)]] <- c(
        list(raw = raw_data),
        likert_data
      )
      
      if (verbose) cat("Simulated: N =", N, "| T =", T, "\n")
    }
  }
  
  data_list
}
