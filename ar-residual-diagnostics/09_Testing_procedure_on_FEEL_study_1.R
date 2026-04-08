
###############################
### MEETING No. 09 ############
###############################

library(dplyr)
library(mgcv)
library(rEDM)
library(nortsTest)
library(ecp)
library(tseries)
library(forecast)
library(lme4)
library(lmerTest)
library(ggplot2)

# library(Hmisc)
# library(vars)
# library(ggplot2)
# library(RColorBrewer)
# library(ggpubr)
# library(tseries)
# library(patchwork)
# library(lme4)
# library(lmerTest)

# set wd 
setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/Research Exchange/PROJECT/Meeting x Pete - Presentations/09")

# get EDM helper
source("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data/SIM_ESM/helpers/HELPER_EDM_tests.R")

# constants
ALPHA_LEVEL <- 0.05 / 18

# load data
data <- read.csv("../05/data_FEEL_Study_1.csv")

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

# Setup progress bar and output data frame
total_iterations <- length(uuids.sub) * length(emotions)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration <- 0
# Create output data frame
results_df <- data.frame(
  UUID = character(),
  emotion = character(),
  p.box = numeric(),
  p.arch.lm = numeric(),
  p.norm.SW = numeric(),
  p.norm.L = numeric(),
  bds.num = numeric(),
  theta_opt = numeric(),
  edf = numeric(),
  cp.num = numeric(),
  stringsAsFactors = FALSE
)

start <- Sys.time()

for (uuid in uuids.sub) {
  temp <- data.sub %>% filter(UUID == uuid)
  
  for (e in emotions) {
    iteration <- iteration + 1
    setTxtProgressBar(pb, iteration)
    
    ts <- na.omit(temp[[e]])
    n <- length(ts)
    time <- 1:n
    
    ts.df <- data.frame(
      time = time,
      ts = ts,
      lag = NA_real_
    )
    
    ts.df$lag <- lag(ts.df$ts, n = 1)
    
    ## 0 - FIT VAR model 
    ar.fit <- arima(ts.df$ts, order = c(1, 0, 0))
    ar.res <- residuals(ar.fit)
    
    # 1.1 - Ljung-Box test
    p.box <- tryCatch(
      Box.test(ar.res, lag = 1, type = "Ljung-Box")$p.value,
      error = function(e) NA
    )
    
    # 1.2 - ARCH LM test
    p.arch.lm <- tryCatch(
      arch.test(ar.res, arch = "Lm")$p.value,
      error = function(e) NA
    )
    
    # 1.3 - Normality 
    p.norm.SW <- tryCatch(
      normal.test(ar.res, normality = c("shapiro"))$p.value,
      error = function(e) NA
    )
    
    p.norm.L <- tryCatch(
      lobato.test(ar.res)$p.value,
      error = function(e) NA
    )
    
    # 1.4 - BDS test (BONUS)
    bds.num <- tryCatch({
      bds <- bds.test(ar.res)
      length(which(bds$p.value < (ALPHA_LEVEL / length(bds$p.value))))
    }, error = function(e) NA)
    
    # 3.1 - Optimal theta
    theta_opt <- tryCatch(
      EDM_tests.ts(ts.df$ts, EDM.include = c("E_opt", "SMap"))[["theta_opt"]],
      error = function(e) NA
    )
    
    # 3.2 - EDF of GAM 
    edf <- tryCatch({
      gam <- gam(ts ~ s(lag, bs = "tp", k = 10), method = "ML", data = ts.df)
      summary(gam)$edf
    }, error = function(e) NA)
    
    # 3.3 - Change Point Analysis
    cp.num <- tryCatch({
      alpha <- 0.05
      R <- ceiling(20 / alpha)
      cp.out <- e.divisive(matrix(ts), R = R, sig.lvl = alpha)
      length(which(cp.out$p.values < alpha))
    }, error = function(e) NA)
    
    # Add row to results data frame
    results_df <- rbind(results_df, data.frame(
      UUID = uuid,
      emotion = e,
      p.box = p.box,
      p.arch.lm = p.arch.lm,
      p.norm.SW = p.norm.SW,
      p.norm.L = p.norm.L,
      bds.num = bds.num,
      theta_opt = theta_opt,
      edf = edf,
      cp.num = cp.num,
      stringsAsFactors = FALSE
    ))
  }
}
close(pb)

end <- Sys.time()

print(end-start)

# View results
head(results_df)



## ---------------------------------
# TRY VAS' APPROACH ----------------

library(forecast)
library(mgcv)

## --- Relative threshold: 2 ---

total_iterations <- length(uuids.sub) * length(emotions)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration <- 0

# Storage for results
results.Vas.rel.th <- data.frame()
# Main loop
for (uuid in uuids.sub) {
  
  temp <- data.sub %>% filter(UUID == uuid)
  
  for (e in emotions) {
    
    tryCatch({
      
      iteration <- iteration + 1
      setTxtProgressBar(pb, iteration)
      
      ts <- na.omit(temp[[e]])
      n <- length(ts)
      time <- 1:n
      
      ts.df <- data.frame(
        time = time,
        ts = ts,
        lag = NA_real_
      )
      
      ts.df$lag <- dplyr::lag(ts.df$ts, n = 1)
      
      # Remove NAs
      ts.df_clean <- na.omit(ts.df)
      n_clean <- nrow(ts.df_clean)
      
      # ===== STEP 2: Fit AR model =====
      ar_model <- tryCatch({
        auto.arima(ts.df_clean$ts,
                   max.p = 3,
                   max.d = 0,
                   max.q = 0,
                   seasonal = FALSE,
                   stepwise = TRUE)
      }, error = function(e) {

        return(NULL)
      })
      
      if (is.null(ar_model)) {
        
        next
      }
      
      ar_fitted <- fitted(ar_model)
      ar_residuals <- ts.df_clean$ts - ar_fitted
      
      # ===== STEP 3: Fit GAM model =====
      gam_model <- tryCatch({
        
        k_value <- min(10, n_clean - 1)
        
        gam(ts ~ s(lag, bs = "tp", k = k_value),
            method = "ML",
            data = ts.df_clean)
        
      }, error = function(e) {
        
        return(NULL)
      })
      
      if (is.null(gam_model)) {
        next
      }
      
      gam_fitted <- fitted(gam_model)
      gam_residuals <- ts.df_clean$ts - gam_fitted
      
      # ===== STEP 4: Flag observations =====
      ar_std_resid <- ar_residuals / sd(ar_residuals)
      gam_std_resid <- gam_residuals / sd(gam_residuals)
      
      ar_flagged <- which(abs(ar_std_resid) > 2)
      gam_flagged <- which(abs(gam_std_resid) > 2)
      
      # ===== STEP 5: Calculate overlap =====
      overlap <- intersect(ar_flagged, gam_flagged)
      
      overlap_pct <- if (max(length(ar_flagged),
                             length(gam_flagged)) > 0) {
        length(overlap) /
          max(length(ar_flagged),
              length(gam_flagged)) * 100
      } else {
        NA
      }
      
      # ===== STEP 6: Error metrics =====
      ar_rmse <- sqrt(mean(ar_residuals^2))
      gam_rmse <- sqrt(mean(gam_residuals^2))
      
      ar_mae <- mean(abs(ar_residuals))
      gam_mae <- mean(abs(gam_residuals))
      
      ar_aic <- ar_model$aic
      gam_aic <- AIC(gam_model)
      
      # ===== STEP 7: Store results =====
      results.Vas.rel.th <- rbind(results.Vas.rel.th, data.frame(
        participant = uuid,
        emotion = e,
        n_obs = n_clean,
        
        ar_rmse = ar_rmse,
        ar_mae = ar_mae,
        ar_aic = ar_aic,
        ar_flagged_n = length(ar_flagged),
        
        gam_rmse = gam_rmse,
        gam_mae = gam_mae,
        gam_aic = gam_aic,
        gam_flagged_n = length(gam_flagged),
        
        overlap_n = length(overlap),
        overlap_pct = overlap_pct,
        
        status = "success"
      ))
      
    }, error = function(e) {
      
      results.Vas.rel.th <<- rbind(results.Vas.rel.th, data.frame(
        participant = uuid,
        emotion = e,
        n_obs = NA,
        ar_rmse = NA,
        ar_mae = NA,
        ar_aic = NA,
        ar_flagged_n = NA,
        gam_rmse = NA,
        gam_mae = NA,
        gam_aic = NA,
        gam_flagged_n = NA,
        overlap_n = NA,
        overlap_pct = NA,
        status = "error"
      ))
    })
  }
}

close(pb)



## --- Absolute threshold: 10 ---

total_iterations <- length(uuids.sub) * length(emotions)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration <- 0

# Storage for results
results.Vas.abs.th <- data.frame()
threshold <- 10

# Main loop
for (uuid in uuids.sub) {
  
  temp <- data.sub %>% filter(UUID == uuid)
  
  for (e in emotions) {
    
    tryCatch({
      
      iteration <- iteration + 1
      setTxtProgressBar(pb, iteration)
      
      ts <- na.omit(temp[[e]])
      n <- length(ts)
      time <- 1:n
      
      ts.df <- data.frame(
        time = time,
        ts = ts,
        lag = NA_real_
      )
      
      ts.df$lag <- dplyr::lag(ts.df$ts, n = 1)
      
      # Remove NAs
      ts.df_clean <- na.omit(ts.df)
      n_clean <- nrow(ts.df_clean)
      
      # ===== STEP 2: Fit AR model =====
      ar_model <- tryCatch({
        auto.arima(ts.df_clean$ts,
                   max.p = 3,
                   max.d = 0,
                   max.q = 0,
                   seasonal = FALSE,
                   stepwise = TRUE)
      }, error = function(e) {
        
        return(NULL)
      })
      
      if (is.null(ar_model)) {
        
        next
      }
      
      ar_fitted <- fitted(ar_model)
      ar_residuals <- ts.df_clean$ts - ar_fitted
      
      # ===== STEP 3: Fit GAM model =====
      gam_model <- tryCatch({
        
        k_value <- min(10, n_clean - 1)
        
        gam(ts ~ s(lag, bs = "tp", k = k_value),
            method = "ML",
            data = ts.df_clean)
        
      }, error = function(e) {
        
        return(NULL)
      })
      
      if (is.null(gam_model)) {
        next
      }
      
      gam_fitted <- fitted(gam_model)
      gam_residuals <- ts.df_clean$ts - gam_fitted
      
      # ===== STEP 4: Flag observations =====
      # ar_std_resid <- ar_residuals / sd(ar_residuals)
      # gam_std_resid <- gam_residuals / sd(gam_residuals)
      # 
      # ar_flagged <- which(abs(ar_std_resid) > 2)
      # gam_flagged <- which(abs(gam_std_resid) > 2)
      
      ar_flagged <- which(abs(ar_residuals) > threshold)
      gam_flagged <- which(abs(gam_residuals) > threshold)
      
      # ===== STEP 5: Calculate overlap =====
      overlap <- intersect(ar_flagged, gam_flagged)
      
      overlap_pct <- if (max(length(ar_flagged),
                             length(gam_flagged)) > 0) {
        length(overlap) /
          max(length(ar_flagged),
              length(gam_flagged)) * 100
      } else {
        NA
      }
      
      # ===== STEP 6: Error metrics =====
      ar_rmse <- sqrt(mean(ar_residuals^2))
      gam_rmse <- sqrt(mean(gam_residuals^2))
      
      ar_mae <- mean(abs(ar_residuals))
      gam_mae <- mean(abs(gam_residuals))
      
      ar_aic <- ar_model$aic
      gam_aic <- AIC(gam_model)
      
      # ===== STEP 7: Store results =====
      results.Vas.abs.th <- rbind(results.Vas.abs.th, data.frame(
        participant = uuid,
        emotion = e,
        n_obs = n_clean,
        
        ar_rmse = ar_rmse,
        ar_mae = ar_mae,
        ar_aic = ar_aic,
        ar_flagged_n = length(ar_flagged),
        
        gam_rmse = gam_rmse,
        gam_mae = gam_mae,
        gam_aic = gam_aic,
        gam_flagged_n = length(gam_flagged),
        
        overlap_n = length(overlap),
        overlap_pct = overlap_pct,
        
        status = "success"
      ))
      
    }, error = function(e) {
      
      results.Vas.abs.th <<- rbind(results.Vas.abs.th, data.frame(
        participant = uuid,
        emotion = e,
        n_obs = NA,
        ar_rmse = NA,
        ar_mae = NA,
        ar_aic = NA,
        ar_flagged_n = NA,
        gam_rmse = NA,
        gam_mae = NA,
        gam_aic = NA,
        gam_flagged_n = NA,
        overlap_n = NA,
        overlap_pct = NA,
        status = "error"
      ))
    })
  }
}

close(pb)
