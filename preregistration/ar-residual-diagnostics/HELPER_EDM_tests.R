EDM_tests <- function(ts_data, EDM.include = c("E_opt", "SMap", "pred_decey")) {

  edm_tests <- c(
    "E_opt",
    "rho_theta0",
    "rho_theta_opt",
    "theta_opt",
    "delta_rho",
    "pred_decay"
  )

  cols <- colnames(ts_data)

  res_person <- matrix(NA, nrow = length(edm_tests), ncol = length(cols),
                       dimnames = list(edm_tests, cols))

  # --- Progress bar setup ---
  pb <- txtProgressBar(min = 0, max = length(cols), style = 3)
  i_col <- 0

  for (col in cols) {
    i_col <- i_col + 1
    setTxtProgressBar(pb, i_col)

    ts <- ts_data[[col]]
    ts_clean <- ts[!is.na(ts)]
    N <- length(ts_clean)
    ts_df <- data.frame(residuals = ts_clean)

    tau <- -1

    # FALLBACK
    # if (is.na(tau)) tau <- 1

    if ("E_opt" %in% EDM.include) {
      mid <- floor(N * 0.7)
      lib <- c(1, mid)
      pred <- c(mid + 1, N)

      # maximum embed dimension to test depends on data length
      Tp <- 1

      e <- tryCatch({
        emb_out <- suppressWarnings(EmbedDimension(
          dataFrame = ts_df,
          lib = lib,
          pred = pred,
          maxE = 15,
          Tp = 1,
          # tau = tau,
          columns = 'residuals',
          target = 'residuals',
          noTime = TRUE,
          showPlot = FALSE
        ))
        emb_out$E[which.max(emb_out$rho)]
      }, error = function(e) NA)

      # print(e)
      if (!is.na(e) && ("SMap" %in% EDM.include)) {
        res_person["E_opt", col] <- e

        thetas <- seq(0, 8, by = 0.5)
        rho_theta <- rep(NA, length(thetas))
        # e <- 1

        for (k in seq_along(thetas)) {
          sm <- tryCatch(
            suppressWarnings(SMap(
              dataFrame = ts_df,
              lib = lib,
              pred = pred,
              columns = "residuals",
              target = "residuals",
              E = e,
              # tau = -1,
              Tp = 1,
              theta = thetas[k],
              noTime = TRUE
            )),
            error = function(e) NULL
          )

          if (!is.null(sm)) {
            sm.obs <- sm$predictions$Observations
            sm.pred <- sm$predictions$Predictions
            valid <- is.finite(sm.obs) & is.finite(sm.pred)

            rho <- cor(sm.obs[valid], sm.pred[valid])
            rho_theta[k] <- rho
          }
        }

        if (any(!is.na(rho_theta))) {
          res_person["rho_theta0", col] <- rho_theta[thetas == 0]
          res_person["rho_theta_opt", col] <- max(rho_theta, na.rm = TRUE)
          res_person["theta_opt", col] <- thetas[which.max(rho_theta)]
          res_person["delta_rho", col] <- res_person["rho_theta_opt", col] - res_person["rho_theta0", col]
        }

        if ("pred_decay" %in% EDM.include) {
          # maximum prediction horizon again depends on pred
          Tp_max <- N - min(pred) - 1
          Tp_max <- max(Tp_max, 1)

          pi <- tryCatch(suppressWarnings(PredictInterval(
            dataFrame = ts_df,
            noTime = TRUE,
            columns = "residuals",
            target = "residuals",
            E = e,
            # tau = tau,
            lib = lib,
            pred = pred,
            maxTp = Tp_max,
            showPlot = FALSE
          )), error = function(e) NULL)

          if (!is.null(pi) && nrow(pi) >= 5) {
            decay <- tryCatch({
              fit <- lm(rho ~ Tp, data = pi[1:5, ])
              coef(fit)[2] * 5
            }, error = function(e) NA)
            res_person["pred_decay", col] <- decay
          }
        }
      }
    }
  }

  close(pb)  # close progress bar

  return(res_person)
}

EDM_tests.ts <- function(ts_data, EDM.include = c("E_opt", "SMap", "pred_decay")) {

  edm_tests <- c(
    "E_opt",
    "rho_theta0",
    "rho_theta_opt",
    "theta_opt",
    "delta_rho",
    "pred_decay"
  )

  # Create results vector instead of matrix
  res_person <- setNames(rep(NA_real_, length(edm_tests)), edm_tests)

  ts <- ts_data
  N <- length(ts)
  tau <- -1

  if ("E_opt" %in% EDM.include) {
    mid <- floor(N * 0.7)
    lib <- c(1, mid)
    pred <- c(mid + 1, N)

    Tp <- 1

    e <- tryCatch({
      emb_out <- suppressWarnings(EmbedDimension(
        dataFrame = as.data.frame(ts),
        lib = lib,
        pred = pred,
        maxE = 15,
        Tp = 1,
        columns = 'ts',
        target = 'ts',
        noTime = TRUE,
        showPlot = FALSE
      ))
      emb_out$E[which.max(emb_out$rho)]
    }, error = function(e) NA)

    if (!is.na(e) && ("SMap" %in% EDM.include)) {
      res_person["E_opt"] <- e

      thetas <- seq(0, 8, by = 0.5)
      rho_theta <- rep(NA, length(thetas))

      for (k in seq_along(thetas)) {
        sm <- tryCatch(
          suppressWarnings(SMap(
            dataFrame = as.data.frame(ts),
            lib = lib,
            pred = pred,
            columns = "ts",
            target = "ts",
            E = e,
            Tp = 1,
            theta = thetas[k],
            noTime = TRUE
          )),
          error = function(e) NULL
        )

        if (!is.null(sm)) {
          sm.obs <- sm$predictions$Observations
          sm.pred <- sm$predictions$Predictions
          valid <- is.finite(sm.obs) & is.finite(sm.pred)

          rho <- cor(sm.obs[valid], sm.pred[valid])
          rho_theta[k] <- rho
        }
      }

      if (any(!is.na(rho_theta))) {
        res_person["rho_theta0"] <- rho_theta[thetas == 0]
        res_person["rho_theta_opt"] <- max(rho_theta, na.rm = TRUE)
        res_person["theta_opt"] <- thetas[which.max(rho_theta)]
        res_person["delta_rho"] <- res_person["rho_theta_opt"] - res_person["rho_theta0"]
      }

      if ("pred_decay" %in% EDM.include) {
        Tp_max <- N - min(pred) - 1
        Tp_max <- max(Tp_max, 1)

        pi <- tryCatch(suppressWarnings(PredictInterval(
          dataFrame = as.data.frame(ts),
          noTime = TRUE,
          columns = "ts",
          target = "ts",
          E = e,
          lib = lib,
          pred = pred,
          maxTp = Tp_max,
          showPlot = FALSE
        )), error = function(e) NULL)

        if (!is.null(pi) && nrow(pi) >= 5) {
          decay <- tryCatch({
            fit <- lm(rho ~ Tp, data = pi[1:5, ])
            coef(fit)[2] * 5
          }, error = function(e) NA)
          res_person["pred_decay"] <- decay
        }
      }
    }
  }

  return(res_person)
}

