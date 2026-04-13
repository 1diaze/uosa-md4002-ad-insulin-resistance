#!/usr/bin/env Rscript --vanilla

run_loo_lmer <- function(
  input_csv,
  output_csv,
  formula_text,
  term_of_interest,
  id_var = "RID",
  base_dir = file.path("/sharedscratch", Sys.getenv("USER"), "md4002-elena-adir"),
  n_cores = NULL,
  use_parallel = TRUE
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(lme4)
    library(parallel)
    library(performance)
  })

  setwd(base_dir)

  df <- read.csv(input_csv, check.names = FALSE)

  bad_names <- which(is.na(names(df)) | trimws(names(df)) == "")
  if (length(bad_names) > 0) {
    names(df)[bad_names] <- paste0("unnamed_", seq_along(bad_names))
  }

  if (!id_var %in% names(df)) {
    stop(sprintf("ID variable '%s' not found in data.", id_var))
  }

  if (is.na(formula_text) || !nzchar(trimws(formula_text))) {
    stop("formula_text is empty.")
  }

  model_formula <- as.formula(formula_text)

  vars_in_formula <- all.vars(model_formula)
  missing_vars <- setdiff(vars_in_formula, names(df))
  if (length(missing_vars) > 0) {
    stop(
      sprintf(
        "Variables missing from data: %s",
        paste(missing_vars, collapse = ", ")
      )
    )
  }

  df[[id_var]] <- as.factor(df[[id_var]])
  subject_levels <- levels(df[[id_var]])
  subject_rows <- split(seq_len(nrow(df)), df[[id_var]])

  ctrl_fast <- lmerControl(
    optimizer = "nloptwrap",
    calc.derivs = FALSE
  )

  fit_full <- lmer(
    formula = model_formula,
    data = df,
    REML = FALSE,
    control = ctrl_fast
  )

  fixef_names <- names(fixef(fit_full))
  if (!term_of_interest %in% fixef_names) {
    stop(
      sprintf(
        "Term '%s' not found in fixed effects.\nAvailable terms: %s",
        term_of_interest,
        paste(fixef_names, collapse = ", ")
      )
    )
  }

  beta_full   <- unname(fixef(fit_full)[term_of_interest])
  se_full     <- unname(sqrt(diag(vcov(fit_full)))[term_of_interest])
  logLik_full <- as.numeric(logLik(fit_full))
  AIC_full    <- AIC(fit_full)
  BIC_full    <- BIC(fit_full)
  sigma_full  <- sigma(fit_full)

  r2_full <- tryCatch(performance::r2(fit_full), error = function(e) NULL)
  r2m_full <- if (is.null(r2_full)) NA_real_ else unname(r2_full$R2_marginal)
  r2c_full <- if (is.null(r2_full)) NA_real_ else unname(r2_full$R2_conditional)

  icc_full <- tryCatch(
    unname(performance::icc(fit_full)$ICC_adjusted),
    error = function(e) NA_real_
  )

  rmse_full <- tryCatch(
    {
      x <- performance::rmse(fit_full)
      if ("RMSE" %in% names(x)) unname(x$RMSE) else as.numeric(x)[1]
    },
    error = function(e) NA_real_
  )

  n_fixed <- length(fixef(fit_full))

  dfbetas_threshold <- 2 / sqrt(length(subject_levels))
  likelihood_threshold <- qchisq(0.95, df = n_fixed)
  aic_threshold <- -2

  if (is.null(n_cores)) {
    n_cores <- max(
      1L,
      as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")) - 1L
    )
  }

  loo_fun <- function(s) {

    drop_rows <- subject_rows[[s]]
    keep_rows <- setdiff(seq_len(nrow(df)), drop_rows)
    df_sub <- df[keep_rows, , drop = FALSE]

    fit_loo <- tryCatch(
      lmer(
        formula = model_formula,
        data = df_sub,
        REML = FALSE,
        control = lmerControl(
          optimizer = "nloptwrap",
          calc.derivs = FALSE
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit_loo)) {
      return(tibble(
        subject = as.character(s),
        n_rows_removed = length(drop_rows),

        beta_full = beta_full,
        se_full = se_full,
        beta_loo = NA_real_,
        delta_beta = NA_real_,
        pct_change_beta = NA_real_,
        dfbetas_term = NA_real_,

        logLik_full = logLik_full,
        logLik_loo = NA_real_,
        delta_logLik = NA_real_,
        likelihood_distance = NA_real_,

        AIC_full = AIC_full,
        AIC_loo = NA_real_,
        delta_AIC = NA_real_,

        BIC_full = BIC_full,
        BIC_loo = NA_real_,
        delta_BIC = NA_real_,

        sigma_full = sigma_full,
        sigma_loo = NA_real_,
        delta_sigma = NA_real_,

        R2_marginal_full = r2m_full,
        R2_marginal_loo = NA_real_,
        delta_R2_marginal = NA_real_,

        R2_conditional_full = r2c_full,
        R2_conditional_loo = NA_real_,
        delta_R2_conditional = NA_real_,

        ICC_full = icc_full,
        ICC_loo = NA_real_,
        delta_ICC = NA_real_,

        RMSE_full = rmse_full,
        RMSE_loo = NA_real_,
        delta_RMSE = NA_real_,

        influential_dfbetas = NA,
        influential_likelihood = NA,
        influential_aic = NA,

        fit_failed = TRUE
      ))
    }

    beta_loo   <- unname(fixef(fit_loo)[term_of_interest])
    logLik_loo <- as.numeric(logLik(fit_loo))
    AIC_loo    <- AIC(fit_loo)
    BIC_loo    <- BIC(fit_loo)
    sigma_loo  <- sigma(fit_loo)

    r2_loo <- tryCatch(performance::r2(fit_loo), error = function(e) NULL)
    r2m_loo <- if (is.null(r2_loo)) NA_real_ else unname(r2_loo$R2_marginal)
    r2c_loo <- if (is.null(r2_loo)) NA_real_ else unname(r2_loo$R2_conditional)

    icc_loo <- tryCatch(
      unname(performance::icc(fit_loo)$ICC_adjusted),
      error = function(e) NA_real_
    )

    rmse_loo <- tryCatch(
      {
        x <- performance::rmse(fit_loo)
        if ("RMSE" %in% names(x)) unname(x$RMSE) else as.numeric(x)[1]
      },
      error = function(e) NA_real_
    )

    dfbetas_term <- (beta_full - beta_loo) / se_full
    likelihood_distance <- 2 * (logLik_full - logLik_loo)
    delta_AIC <- AIC_loo - AIC_full

    tibble(
      subject = as.character(s),
      n_rows_removed = length(drop_rows),

      beta_full = beta_full,
      se_full = se_full,
      beta_loo = beta_loo,
      delta_beta = beta_full - beta_loo,
      pct_change_beta = ifelse(beta_full == 0, NA_real_, 100 * (beta_loo - beta_full) / beta_full),
      dfbetas_term = dfbetas_term,

      logLik_full = logLik_full,
      logLik_loo = logLik_loo,
      delta_logLik = logLik_full - logLik_loo,
      likelihood_distance = likelihood_distance,

      AIC_full = AIC_full,
      AIC_loo = AIC_loo,
      delta_AIC = delta_AIC,

      BIC_full = BIC_full,
      BIC_loo = BIC_loo,
      delta_BIC = BIC_loo - BIC_full,

      sigma_full = sigma_full,
      sigma_loo = sigma_loo,
      delta_sigma = sigma_loo - sigma_full,

      R2_marginal_full = r2m_full,
      R2_marginal_loo = r2m_loo,
      delta_R2_marginal = r2m_loo - r2m_full,

      R2_conditional_full = r2c_full,
      R2_conditional_loo = r2c_loo,
      delta_R2_conditional = r2c_loo - r2c_full,

      ICC_full = icc_full,
      ICC_loo = icc_loo,
      delta_ICC = icc_loo - icc_full,

      RMSE_full = rmse_full,
      RMSE_loo = rmse_loo,
      delta_RMSE = rmse_loo - rmse_full,

      influential_dfbetas = abs(dfbetas_term) > dfbetas_threshold,
      influential_likelihood = likelihood_distance > likelihood_threshold,
      influential_aic = delta_AIC < aic_threshold,

      fit_failed = FALSE
    )
  }

  if (use_parallel) {
    loo_results_list <- mclapply(subject_levels, loo_fun, mc.cores = n_cores)
  } else {
    loo_results_list <- lapply(subject_levels, loo_fun)
  }

  loo_results <- bind_rows(loo_results_list)

  full_model_metrics <- tibble(
    formula = formula_text,
    term_of_interest = term_of_interest,
    id_var = id_var,
    beta_full = beta_full,
    se_full = se_full,
    n_obs = nobs(fit_full),
    n_subjects = length(subject_levels),
    AIC = AIC_full,
    BIC = BIC_full,
    logLik = logLik_full,
    sigma = sigma_full,
    R2_marginal = r2m_full,
    R2_conditional = r2c_full,
    ICC = icc_full,
    RMSE = rmse_full,
    dfbetas_threshold = dfbetas_threshold,
    likelihood_threshold = likelihood_threshold,
    aic_threshold = aic_threshold
  )

  write.csv(loo_results, output_csv, row.names = FALSE)

  metrics_output_csv <- sub("\\.csv$", "_full_metrics.csv", output_csv)
  if (identical(metrics_output_csv, output_csv)) {
    metrics_output_csv <- paste0(output_csv, "_full_metrics.csv")
  }
  write.csv(full_model_metrics, metrics_output_csv, row.names = FALSE)

  return(list(
    loo_results = loo_results,
    full_model_metrics = full_model_metrics,
    fit_full = fit_full
  ))
}