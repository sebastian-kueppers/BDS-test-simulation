
#####################################
### AR RESIDUAL DIAGNOSTICS #########
#####################################

library(dplyr)
library(FKF)
library(tseries) 
library(tidyverse)
library(future)
library(furrr)
library(FinTS)
library(ecp)
library(effectsize)
library(FKF)
library(foreign)
library(rjags)
library(coda)
library(posterior)


# set wd 
# setwd("...")
setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/Research Exchange/PROJECT/Preregistration/github/ar-residual-diagnostics")

# load data
data <- read.csv("data_FEEL_Study_1.csv")

# get helper function
source("get_ARWN_residuals.R")
source("HELPER_EDM_tests.R")


# get all uuids
uuids <- unique(data$UUID) # 179 participants

# draw 1/3 of sample for testings
set.seed(2310)
uuids.sub <- sample(uuids, round(length(uuids) / 3)) # 60

# check if seeding worked
uuids.sub[3] # "5b4774af-4900-4d02-80ad-dc3a334e868f"
uuids.sub[23] # "4ea6afa1-e202-45eb-80c5-ad90f1ce6ce5"

# retrieve sub data
data.sub <- data[data$UUID %in% uuids.sub, ]

# get affect items
items.neg <- c("ANG_ES", "SAD_ES", "STR_ES")
items.pos <- c("CONF_ES", "HAP_ES", "RLX_ES") 
emotions <- c(items.pos, items.neg) 

## --------------------------------
# PROCEDURE -----------------------
## --------------------------------

## --------------------------------
# Setup parallel ------------------
# 
n_cores <- availableCores() - 1
plan(multisession, workers = n_cores)

# Prepare data for parallel processing

task_data <- expand_grid(
  uuid = uuids.sub,
  emotion = emotions
) %>%
  mutate(
    # ts = map2(uuid, emotion, function(u, e) {
    #   temp <- data.sub %>% filter(UUID == u)
    #   na.omit(temp[[e]])
    # })
    
    # with overnight lag:
    ts = map2(uuid, emotion, function(u, e) {
      temp <- data.sub %>% 
        filter(UUID == u) %>%
        arrange(Date_Local, Time_Local) 
      
      y <- temp[[e]]
      dates <- temp$Date_Local
      
      ts_with_overnight <- c()
      for (i in seq_along(y)) {
        if (i > 1 && dates[i] != dates[i - 1]) {
          ts_with_overnight <- c(ts_with_overnight, NA)  # Overnight NA
        }
        ts_with_overnight <- c(ts_with_overnight, y[i])
      }
      
      ts_with_overnight  # NAs from missing beeps remain
    })
  )

# Some descriptives
# n_mean <- mean(map_dbl(task_data$ts, length)) # 163.5139

## -------------------------------
# Helper function ----------------

fit_ar1wn_residuals_parallel <- function(ts, uuid, emotion) {
  tryCatch({
    
    # Check inclusion criteria (flag only, do not skip analysis)
    ts_complete <- ts[!is.na(ts)]
    excluded <- FALSE
    exclusion_reason <- NA
    
    if (length(ts_complete) < 50) {
      excluded <- TRUE
      exclusion_reason <- "fewer than 50 complete observations"
    } else if (sd(ts_complete) == 0) {
      excluded <- TRUE
      exclusion_reason <- "SD = 0"
    } else if (length(unique(ts_complete)) < 3) {
      excluded <- TRUE
      exclusion_reason <- "fewer than 3 response categories"
    }
    
    freq_result  <- NULL
    bayes_result <- NULL
    
    if (!excluded) {
      # Frequentist
      freq_result <- fit_ar1wn_residuals(ts, verbose = FALSE)
      # Bayesian
      bayes_result <- fit_ar1wn_residuals_bayesian(ts, verbose = TRUE)
    }

    # Return both
    list(
      uuid = uuid,
      emotion = emotion,
      frequentist = freq_result,
      bayesian = bayes_result,
      excluded = excluded,
      excusion_reason = exclusion_reason,
      success = TRUE
    )
  }, error = function(e) {
    list(
      uuid = uuid,
      emotion = emotion,
      error = e$message,
      success = FALSE
    )
  })
}

## ------------------------------
# Run parallel loop -------------

start_time <- Sys.time()

# Single
results_parallel <- task_data %>%
  mutate(
    result = future_pmap(
      list(ts, uuid, emotion),
      fit_ar1wn_residuals_parallel,
      .progress = TRUE
    )
  )

# test_results <- task_data %>%
#   slice(1:3) %>%
#   mutate(
#     result = future_pmap(
#       list(ts, uuid, emotion),
#       fit_ar1wn_residuals_parallel,
#       .progress = TRUE
#     )
#   )

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

log_message <- sprintf(
  "[%s] Dauer: %.2f Minuten\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  as.numeric(duration)
)
cat(log_message, file = "log.txt", append = TRUE)

# Save
saveRDS(results_parallel, file = "ARWN_results_19-05.rds")

plan(sequential)


# Read data
results <- readRDS("ARWN_results.rds")

## --------------------------------
# CONVERGENCE AND HEYWOOD CASES ---
## --------------------------------

heywood_ivar <- c()
heywood_evar <- c()
bayes_convergence_failed <- c()

for (i in 1:length(results$result)) {
  freq <- results$result[[i]]$frequentist
  bayesian <- results$result[[i]]$bayesian
  if (!is.null(bayesian)) {
    bayes_convergence_failed <- c(bayes_convergence_failed, 
                                  bayesian$convergence_failed)
  }
  
  heywood_ivar <- c(heywood_ivar, freq$heywood_ivar)
  heywood_evar <- c(heywood_evar, freq$heywood_evar)
}

sum(bayes_convergence_failed) # 360

# Bayesian estimations did not converge at all. Problem: estimations of ivar and
# evar are highly correlated, with up to r = -.91.

sum(heywood_ivar | heywood_evar) # 181

# 181/360 = 50% of frequentist estimations show Heywood cases ... 



## --------------------------------
# ANALYSIS ------------------------
## --------------------------------

source("bootstrap_BDS.R")

# Excluding freq. Heywood cases likely introduces a selection bias because time
# series with lower AR effect are known to produce more such cases (Schuurman et al., 2015)

# Hence, we will continue with Bayesian cases that did not fail in estimation.

# Filter to only successful Bayesian fits
valid_results <- results %>%
  mutate(
    has_bayesian = map_lgl(result, ~!is.null(.$bayesian))
  ) %>%
  filter(has_bayesian)

# Initialize empty list
BDS_results <- list()

# Run BDS tests sequentially with progress bar
n_total <- nrow(valid_results)

for (i in 1:n_total) {
  uuid <- valid_results$uuid[i]
  emotion <- valid_results$emotion[i]
  ts <- valid_results$ts[i]
  residuals <- valid_results$result[[i]]$bayesian$residuals
  
  bds_result <- bds_test_bootstrap(residuals)
  
  # control: BDS result without controlling for measurement error
  ar.fit <- arima(ts[[1]], order = c(1, 0, 0))
  ar.res <- residuals(ar.fit)
  bds_result_w_me <- bds_test_bootstrap(ar.res)
  
  BDS_results[[i]] <- list(
    uuid = uuid,
    emotion = emotion,
    ts = ts,
    bds_result = bds_result,
    bds_result_w_me = bds_result_w_me
  )
  
  # Print progress
  cat("\r", sprintf("Progress: %d/%d (%.1f%%)", i, n_total, i/n_total*100), sep = "")
  flush.console()
}

cat("\n")  # New line after progress bar

# Convert to tibble
BDS_results <- map_df(BDS_results, ~tibble(
  uuid = .$uuid,
  emotion = .$emotion,
  bds_result = list(.$bds_result),
  bds_result_w_me = list(.$bds_result_w_me)
))

# Proportion of significant BDS tests
prop_sig <- mean(map_dbl(BDS_results$bds_result, ~.$bds_pvalue_empirical) < 0.05) #0.5444
prop_sig_w_me <- mean(map_dbl(BDS_results$bds_result_w_me, ~.$bds_pvalue_empirical) < 0.05) #0.503


# Check per participant
# Extract p-values and add to BDS_results
BDS_results_with_pval <- BDS_results %>%
  mutate(
    bds_pvalue_empirical = map_dbl(bds_result, ~.$bds_pvalue_empirical),
    significant = bds_pvalue_empirical < 0.05,
    bds_pvalue_empirical_w_me = map_dbl(bds_result_w_me, ~.$bds_pvalue_empirical),
    significant_w_me = bds_pvalue_empirical_w_me < 0.05
  )

# Does controlling for measurement error have an effect on BDS results?
table(BDS_results_with_pval$significant, 
      BDS_results_with_pval$significant_w_me)

#       FALSE TRUE
# FALSE   145   17
# TRUE     35  163

# ========== Summary by UUID ==========
bds_by_uuid <- BDS_results_with_pval %>%
  group_by(uuid) %>%
  summarise(
    n_emotions = n(),
    n_significant = sum(significant),
    prop_significant = n_significant / n_emotions,
    .groups = "drop"
  ) %>%
  arrange(desc(prop_significant))

mean(bds_by_uuid$n_significant) # 3.266667
sd(bds_by_uuid$n_significant) # 1.725842


## --------------------------------
# VALIDATION AND COMPLEXITY -------
## --------------------------------

results_diag <- list()
results_complex <- list()

# Loop over results
for (i in seq_len(nrow(results))) {
  print(i)
  
  # Extract key information
  uuid <- results$uuid[i]
  emotion <- results$emotion[i]
  ts <- results$ts[i]
  
  # Extract residuals
  residuals_ts <- results$result[[i]]$bayesian$residuals
  
  # Extract BDS result
  bds_match <- BDS_results_with_pval %>%
    filter(uuid == !!uuid & emotion == !!emotion)
  
  BDS.p <- if (nrow(bds_match) > 0) {
    bds_match$bds_pvalue_empirical[1]
  } else {
    NA
  }
  
  # Skip if residuals are missing or insufficient
  if (is.null(residuals_ts) || length(residuals_ts) < 10) {
    next
  }
  
  # Create time index
  time_idx <- 1:length(residuals_ts)
  
  # ========== GAM MODELLING ==========
  
  # Create lagged residuals
  residuals_df <- data.frame(
    residuals = residuals_ts,
    residuals_lag = c(NA, residuals_ts[-length(residuals_ts)]),
    time = time_idx
  )
  
  # Remove NA rows
  residuals_df <- residuals_df %>% filter(!is.na(residuals_lag))
  
  edf.lag <- NA
  edf.t <- NA
  
  tryCatch({
    gam.lag <- gam(residuals ~ s(residuals_lag, bs = "tp", k = 10), 
                   method = "ML", data = residuals_df)
    edf.lag <- summary(gam.lag)$edf
  }, error = function(e) {
    message(paste("GAM lag error for", uuid, ":", e$message))
  })
  
  tryCatch({
    gam.t <- gam(residuals ~ s(time, bs = "tp", k = 10), 
                 method = "ML", data = residuals_df)
    edf.t <- summary(gam.t)$edf
  }, error = function(e) {
    message(paste("GAM time error for", uuid, ":", e$message))
  })
  
  # ========== S-MAP ANALYSIS (EDM) ==========
  
  E_opt <- NA
  theta_opt <- NA
  
  tryCatch({
    edm <- EDM_tests(data.frame(residuals = residuals_ts), 
                     EDM.include = c("E_opt", "SMap"))
    E_opt <- edm["E_opt", "residuals"]
    theta_opt <- edm["theta_opt", "residuals"]
  }, error = function(e) {
    message(paste("EDM error for", uuid, ":", e$message))
  })
  
  # ========== LJUNG-BOX TEST ==========
  
  ljung_box.p <- NA
  
  tryCatch({
    ljung_box <- Box.test(residuals_ts, lag = 10, type = "Ljung-Box")
    ljung_box.p <- ljung_box$p.value
  }, error = function(e) {
    message(paste("Ljung-Box error for", uuid, ":", e$message))
  })
  
  # ========== ARCH LM TEST ==========
  
  arch_lm.p <- NA
  
  tryCatch({
    arch_lm <- ArchTest(residuals_ts, lags = 10)
    arch_lm.p <- arch_lm$p.value
  }, error = function(e) {
    message(paste("ARCH LM error for", uuid, ":", e$message))
  })
  
  # ========== SHAPIRO-WILK TEST ==========
  
  shapiro_wilk.p <- NA
  
  tryCatch({
    shapiro_wilk <- shapiro.test(residuals_ts)
    shapiro_wilk.p <- shapiro_wilk$p.value
  }, error = function(e) {
    message(paste("Shapiro-Wilk error for", uuid, ":", e$message))
  })
  
  # ========== CHANGE POINT ANALYSIS ==========
  
  changepoint.n <- NA
  
  tryCatch({
    alpha <- 0.05
    R <- ceiling(20 / alpha)
    cp.out <- e.divisive(matrix(residuals_ts), R = R, sig.lvl = alpha)
    changepoint.n <- length(which(cp.out$p.values < alpha))
  }, error = function(e) NA)
  
  # ========== STORE RESULTS ==========
  
  # Create diagnostic test row
  results_diag[[i]] <- data.frame(
    uuid = uuid,
    emotion = emotion,
    # ts = c(ts),
    n = length(residuals_ts),
    ljung_box.p = ljung_box.p,
    arch_lm.p = arch_lm.p,
    shapiro_wilk.p = shapiro_wilk.p,
    BDS.p = BDS.p,
    stringsAsFactors = FALSE
  )
  
  row.names(diag_row) <- NULL
  
  # Create complexity measures row
  results_complex[[i]] <- data.frame(
    uuid = uuid,
    emotion = emotion,
    # ts = ts,
    n = length(residuals_ts),
    edf.lag = edf.lag,
    edf.t = edf.t,
    E_opt = E_opt,
    theta_opt = theta_opt,
    changepoint.n = changepoint.n,
    BDS.p = BDS.p,
    stringsAsFactors = FALSE
  )
}



# After the loop, bind all rows at once
results_diag <- do.call(rbind, results_diag)
results_complex <- do.call(rbind, results_complex)

# Reset row names
rownames(results_diag) <- NULL
rownames(results_complex) <- NULL

# Save dataframes
saveRDS(results_diag, file = "results_diag.rds")
saveRDS(results_complex, file = "results_complex.rds")


## --------------------------------
# CHECK RESULTS -------------------
## --------------------------------

## --------------------------------
# Validating BDS Test -------------

# Conduct Chi square tests to test whether rows with signficant BDS test are more 
# likely to have significant Ljung-Box, ARCH LM, and Shapiro-Wilk tests.

# Create binary indicators for significance (p < 0.05)
results_diag <- results_diag %>%
  mutate(
    BDS_sig = BDS.p < 0.05,
    ljung_box_sig = ljung_box.p < 0.05,
    arch_lm_sig = arch_lm.p < 0.05,
    shapiro_wilk_sig = shapiro_wilk.p < 0.05
  )

# Ljung-Box vs BDS
ct_ljung <- table(results_diag$BDS_sig, results_diag$ljung_box_sig)
chi_ljung <- chisq.test(ct_ljung)
chi_ljung # significant

# ARCH LM vs BDS
ct_arch <- table(results_diag$BDS_sig, results_diag$arch_lm_sig)
chi_arch <- chisq.test(ct_arch)
chi_arch # significant

# Shapiro-Wilk vs BDS
ct_shapiro <- table(results_diag$BDS_sig, results_diag$shapiro_wilk_sig)
chi_shapiro <- chisq.test(ct_shapiro)
chi_shapiro # non-significant


## --------------------------------
# Complexity Markers --------------

# Conduct Mann-Whitney U tests to test whether rows that are significant on the 
# BDS test differ in their values from rows that are non-significant

# Create binary indicator for BDS significance
results_complex <- results_complex %>%
  mutate(
    BDS_sig = BDS.p < 0.05
  )

summary_stats <- results_complex %>%
  group_by(BDS_sig) %>%
  summarise(
    n = n(),
    # EDF Lag
    edf.lag_mean = mean(edf.lag, na.rm = TRUE),
    edf.lag_sd = sd(edf.lag, na.rm = TRUE),
    edf.lag_median = median(edf.lag, na.rm = TRUE),
    # EDF Time
    edf.t_mean = mean(edf.t, na.rm = TRUE),
    edf.t_sd = sd(edf.t, na.rm = TRUE),
    edf.t_median = median(edf.t, na.rm = TRUE),
    # Theta Optimal
    theta_opt_mean = mean(theta_opt, na.rm = TRUE),
    theta_opt_sd = sd(theta_opt, na.rm = TRUE),
    theta_opt_median = median(theta_opt, na.rm = TRUE),
    # Changepoint Count
    changepoint.n_mean = mean(changepoint.n, na.rm = TRUE),
    changepoint.n_sd = sd(changepoint.n, na.rm = TRUE),
    changepoint.n_median = median(changepoint.n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    BDS_sig = ifelse(BDS_sig, "BDS p < 0.05", "BDS p ≥ 0.05")
  )


# ========== MANN-WHITNEY U TESTS ==========

mw_edf_lag <- wilcox.test(
  results_complex$edf.lag[results_complex$BDS_sig == TRUE],
  results_complex$edf.lag[results_complex$BDS_sig == FALSE],
  alternative = "two.sided"
)
mw_edf_lag # significant

rank_biserial(
  results_complex$edf.lag[results_complex$BDS_sig == TRUE],
  results_complex$edf.lag[results_complex$BDS_sig == FALSE]
) # r = 0.32


mw_edf_t <- wilcox.test(
  results_complex$edf.t[results_complex$BDS_sig == TRUE],
  results_complex$edf.t[results_complex$BDS_sig == FALSE],
  alternative = "two.sided"
)
mw_edf_t # non-significant

rank_biserial(
  results_complex$edf.t[results_complex$BDS_sig == TRUE],
  results_complex$edf.t[results_complex$BDS_sig == FALSE]
) # r < 0.01


mw_theta_opt <- wilcox.test(
  results_complex$theta_opt[results_complex$BDS_sig == TRUE],
  results_complex$theta_opt[results_complex$BDS_sig == FALSE],
  alternative = "two.sided"
)
mw_theta_opt # non-significant

rank_biserial(
  results_complex$theta_opt[results_complex$BDS_sig == TRUE],
  results_complex$theta_opt[results_complex$BDS_sig == FALSE]
) # r = -0.04


mw_changepoint <- wilcox.test(
  results_complex$changepoint.n[results_complex$BDS_sig == TRUE],
  results_complex$changepoint.n[results_complex$BDS_sig == FALSE],
  alternative = "two.sided"
)
mw_changepoint # significant

rank_biserial(
  results_complex$changepoint.n[results_complex$BDS_sig == TRUE],
  results_complex$changepoint.n[results_complex$BDS_sig == FALSE]
) # r = 0.16

## --------------------------------
# EXPLORE: /WO MEASUREMENT ERROR CORRECTION
## --------------------------------

# Iterate over results

