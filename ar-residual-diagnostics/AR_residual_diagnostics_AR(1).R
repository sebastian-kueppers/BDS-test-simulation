## --------------------------------
# AR RESIDUAL DIAGNOSTICS ---------
## --------------------------------

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
library(lme4)
library(mgcv)
library(rEDM)


# set wd
# setwd("...")
setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/Research Exchange/PROJECT/Preregistration/github/ar-residual-diagnostics")

# Data can be given out on request!
# load data
data <- read.csv("data_FEEL_Study_1.csv")

# get helper function
source("get_ARWN_residuals.R")
source("HELPER_EDM_tests.R")
source("bootstrap_BDS.R")

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
# Setup Parallel ------------------

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

## --------------------------------
# Run BDS Testing Procedure -------

BDS_results <- list()

pb <- txtProgressBar(min = 0, max = nrow(task_data), style = 3)
pb_count <- 0

set.seed(2310)

for (i in 1:nrow(task_data)) {
  uuid <- task_data$uuid[i]
  emotion <- task_data$emotion[i]
  ts <- task_data$ts[i][[1]]
  
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
  
  ar.fit <- arima(ts, order = c(1, 0, 0))
  ar.res <- residuals(ar.fit)
  ar.res.wo.na <- ar.res[!is.na(ar.res)]
  bds_result_w_me <- bds_test_bootstrap(ar.res.wo.na)
  
  BDS_results[[i]] <- list(
    uuid = uuid,
    emotion = emotion,
    ts = ts,
    bds_result_w_me = bds_result_w_me,
    bds_statistic = bds_result_w_me$bds_statistic,
    bds_pvalue_parametric = bds_result_w_me$bds_pvalue_parametric,
    bds_pvalue_empirical = bds_result_w_me$bds_pvalue_empirical,
    residuals = ar.res,
    excluded = excluded,
    exclusion_reason = exclusion_reason
  )
  
  pb_count <- pb_count + 1
  setTxtProgressBar(pb, pb_count)
}


# Save
saveRDS(BDS_results, file = "BDS_results_AR_02-06.rds")

# Read data
BDS_results <- readRDS("BDS_results_AR_02-06.rds")


# Convert to tibble
BDS_results.tibble <- map_df(BDS_results, ~tibble(
  uuid = .$uuid,
  emotion = .$emotion,
  excluded = .$excluded,
  exclusion_reason = .$exclusion_reason,
  bds_result_w_me = list(.$bds_result_w_me)
))

# Proportion of time series failing inclusion criteria
mean(BDS_results.tibble$excluded) # 0 -> all meet criteria

# Proportion of significant BDS tests
mean(map_dbl(BDS_results.tibble$bds_result_w_me, ~.$bds_pvalue_empirical) < 0.05) #0.488129


# Check per participant
# Extract p-values and add to BDS_results
BDS_results_with_pval <- BDS_results.tibble %>%
  mutate(
    bds_pvalue_empirical = map_dbl(bds_result_w_me, ~.$bds_pvalue_empirical),
    significant = bds_pvalue_empirical < 0.05,
  )


## --------------------------------
# Summary by UUID and Valence -----

emotion_valence <- tibble(
  emotion = c(items.pos, items.neg),
  valence  = c(rep("positive", length(items.pos)),
               rep("negative", length(items.neg)))
)

# Join and summary
bds_by_valence <- BDS_results_with_pval %>%
  left_join(emotion_valence, by = "emotion") %>%
  group_by(valence) %>%
  summarise(
    n_emotions      = n(),
    n_significant   = sum(significant),
    prop_significant = mean(significant),
    .groups = "drop"
  )

bds_by_valence

BDS_results_with_pval %>%
  left_join(emotion_valence, by = "emotion") %>%
  glmer(significant ~ valence + (1 | uuid),
        data   = .,
        family = binomial) %>%
  summary()

bds_by_uuid_valence <- BDS_results_with_pval %>%
  left_join(emotion_valence, by = "emotion") %>%
  group_by(uuid, valence) %>%
  summarise(
    n_significant    = sum(significant),
    prop_significant = mean(significant),
    .groups = "drop"
  )

BDS_results_with_pval %>%
  group_by(uuid) %>%
  summarise(n_significant = sum(significant)) %>%
  summarise(
    mean_n_significant = mean(n_significant),
    sd_n_significant   = sd(n_significant)
  )
# mean_n_significant sd_n_significant
# <dbl>            <dbl>
#   1               2.93             1.73

bds_by_uuid_valence %>%
  group_by(valence) %>%
  summarise(
    mean_n_significant   = mean(n_significant),
    sd_n_significant     = sd(n_significant),
    mean_prop_significant = mean(prop_significant),
    sd_prop_significant   = sd(prop_significant)
  )

# valence  mean_n_significant sd_n_significant mean_prop_significant sd_prop_significant
# <chr>                 <dbl>            <dbl>                 <dbl>               <dbl>
#   1 negative               1.55            0.946                 0.517               0.315
# 2 positive               1.38            1.06                  0.461               0.353


## --------------------------------
# VALIDATION AND COMPLEXITY -------
## --------------------------------

results_diag <- list()
results_complex <- list()

# Loop over results
for (i in 1:nrow(task_data)) {
  print(i)
  
  uuid <- task_data$uuid[i]
  emotion <- task_data$emotion[i]
  ts <- task_data$ts[i][[1]]
  
  # Extract residuals
  residuals_ts <- BDS_results[[i]]$residuals
  
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
  
  ## --------------------------------
  # GAM Modelling (Raw TS) ----------
  
  ts_complete <- ts[!is.na(ts)]
  
  ts_df <- data.frame(
    ts      = ts_complete,
    ts_lag  = c(NA, ts_complete[-length(ts_complete)]),
    time    = 1:length(ts_complete)
  ) %>% filter(!is.na(ts_lag))
  
  edf.lag <- NA
  edf.t   <- NA
  
  tryCatch({
    gam.lag <- gam(ts ~ s(ts_lag, bs = "tp", k = 10),
                   method = "ML", data = ts_df)
    edf.lag <- summary(gam.lag)$edf
  }, error = function(e) {
    message(paste("GAM lag error for", uuid, ":", e$message))
  })
  
  tryCatch({
    gam.t <- gam(ts ~ s(time, bs = "tp", k = 10),
                 method = "ML", data = ts_df)
    edf.t <- summary(gam.t)$edf
  }, error = function(e) {
    message(paste("GAM time error for", uuid, ":", e$message))
  })
  
  ## --------------------------------
  # S-Map Analysis (EDM) ------------
  
  E_opt     <- NA
  theta_opt <- NA
  
  tryCatch({
    edm <- EDM_tests(data.frame(residuals = ts),
                     EDM.include = c("E_opt", "SMap"))
    E_opt     <- edm["E_opt",     "residuals"]
    theta_opt <- edm["theta_opt", "residuals"]
  }, error = function(e) {
    message(paste("EDM error for", uuid, ":", e$message))
  })
  
  ## --------------------------------
  # Ljung-Box Test ------------------
  
  ljung_box.p <- NA
  
  tryCatch({
    ljung_box   <- Box.test(residuals_ts, lag = 10, type = "Ljung-Box")
    ljung_box.p <- ljung_box$p.value
  }, error = function(e) {
    message(paste("Ljung-Box error for", uuid, ":", e$message))
  })
  
  ## --------------------------------
  # ARCH LM Test --------------------
  
  arch_lm.p <- NA
  
  tryCatch({
    residuals_clean <- residuals_ts[!is.na(residuals_ts)]
    n_lags <- 10
    if (length(residuals_clean) > n_lags + 5) {
      arch_lm   <- ArchTest(residuals_clean, lags = n_lags)
      arch_lm.p <- arch_lm$p.value
    }
  }, error = function(e) {
    message(paste("ARCH LM error for", uuid, ":", e$message))
  })
  
  ## --------------------------------
  # Shapiro-Wilk Test ---------------
  
  shapiro_wilk.p <- NA
  
  tryCatch({
    shapiro_wilk   <- shapiro.test(residuals_ts)
    shapiro_wilk.p <- shapiro_wilk$p.value
  }, error = function(e) {
    message(paste("Shapiro-Wilk error for", uuid, ":", e$message))
  })
  
  ## --------------------------------
  # Change Point Analysis (Raw TS) --
  
  changepoint.n <- NA
  
  tryCatch({
    alpha         <- 0.05
    R             <- ceiling(20 / alpha)
    cp.out        <- e.divisive(matrix(ts_complete), R = R, sig.lvl = alpha)
    changepoint.n <- length(which(cp.out$p.values < alpha))
  }, error = function(e) NA)
  
  ## --------------------------------
  # Store Results -------------------
  
  results_diag[[i]] <- data.frame(
    uuid           = uuid,
    emotion        = emotion,
    n              = length(residuals_ts),
    ljung_box.p    = ljung_box.p,
    arch_lm.p      = arch_lm.p,
    shapiro_wilk.p = shapiro_wilk.p,
    BDS.p          = BDS.p,
    stringsAsFactors = FALSE
  )
  
  results_complex[[i]] <- data.frame(
    uuid          = uuid,
    emotion       = emotion,
    n             = length(residuals_ts),
    edf.lag       = edf.lag,
    edf.t         = edf.t,
    E_opt         = E_opt,
    theta_opt     = theta_opt,
    changepoint.n = changepoint.n,
    BDS.p         = BDS.p,
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
saveRDS(results_diag, file = "results_diag_02-06.rds")
saveRDS(results_complex, file = "results_complex_02-06.rds")


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


## --------------------------------
# Mann-Whitney U Tests ------------

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
mw_edf_t # significant

rank_biserial(
  results_complex$edf.t[results_complex$BDS_sig == TRUE],
  results_complex$edf.t[results_complex$BDS_sig == FALSE]
) # r < 0.31


mw_theta_opt <- wilcox.test(
  results_complex$theta_opt[results_complex$BDS_sig == TRUE],
  results_complex$theta_opt[results_complex$BDS_sig == FALSE],
  alternative = "two.sided"
)
mw_theta_opt # non-significant

rank_biserial(
  results_complex$theta_opt[results_complex$BDS_sig == TRUE],
  results_complex$theta_opt[results_complex$BDS_sig == FALSE]
) # r = 0.04


mw_changepoint <- wilcox.test(
  results_complex$changepoint.n[results_complex$BDS_sig == TRUE],
  results_complex$changepoint.n[results_complex$BDS_sig == FALSE],
  alternative = "two.sided"
)
mw_changepoint # significant

rank_biserial(
  results_complex$changepoint.n[results_complex$BDS_sig == TRUE],
  results_complex$changepoint.n[results_complex$BDS_sig == FALSE]
) # r = 0.30