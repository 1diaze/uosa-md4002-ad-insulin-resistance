suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

plot_full_model_metrics <- function(
  metrics_all,
  metrics = c("R2_marginal", "R2_conditional", "AIC", "BIC", "ICC", "RMSE", "sigma"),
  outcome_order = NULL
) {
  dat <- metrics_all %>%
    select(any_of(c("outcome", metrics))) %>%
    pivot_longer(
      cols = -outcome,
      names_to = "metric",
      values_to = "value"
    )

  if (!is.null(outcome_order)) {
    dat$outcome <- factor(dat$outcome, levels = outcome_order)
  }

  ggplot(dat, aes(x = outcome, y = value)) +
    geom_point(size = 3) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw() +
    labs(
      x = NULL,
      y = "Metric value",
      title = "Full-model performance across outcomes"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_loo_metric_distributions <- function(
  loo_all,
  metrics = c("dfbetas_term", "delta_AIC", "delta_BIC", "delta_R2_marginal", "delta_R2_conditional", "delta_RMSE"),
  outcome_order = NULL,
  add_jitter = TRUE
) {
  dat <- loo_all %>%
    select(any_of(c("outcome", "subject", metrics))) %>%
    pivot_longer(
      cols = -c(outcome, subject),
      names_to = "metric",
      values_to = "value"
    )

  if (!is.null(outcome_order)) {
    dat$outcome <- factor(dat$outcome, levels = outcome_order)
  }

  p <- ggplot(dat, aes(x = outcome, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw() +
    labs(
      x = NULL,
      y = "LOO change",
      title = "Leave-one-subject-out sensitivity"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if (add_jitter) {
    p <- p + geom_jitter(width = 0.15, alpha = 0.30, size = 0.9)
  }

  p
}

summarise_loo_robustness <- function(loo_all) {
  loo_all %>%
    group_by(outcome) %>%
    summarise(
      median_abs_dfbetas = median(abs(dfbetas_term), na.rm = TRUE),
      max_abs_dfbetas = max(abs(dfbetas_term), na.rm = TRUE),
      mean_abs_dfbetas = mean(abs(dfbetas_term), na.rm = TRUE),
      sd_delta_AIC = sd(delta_AIC, na.rm = TRUE),
      sd_delta_BIC = sd(delta_BIC, na.rm = TRUE),
      sd_delta_R2_marginal = sd(delta_R2_marginal, na.rm = TRUE),
      sd_delta_R2_conditional = sd(delta_R2_conditional, na.rm = TRUE),
      sd_delta_RMSE = sd(delta_RMSE, na.rm = TRUE),
      n_influential_dfbetas = sum(influential_dfbetas, na.rm = TRUE),
      n_failed = sum(fit_failed, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_loo_robustness_summary <- function(
  loo_summary,
  metrics = c("median_abs_dfbetas", "max_abs_dfbetas", "sd_delta_AIC", "sd_delta_R2_marginal"),
  outcome_order = NULL
) {
  dat <- loo_summary %>%
    select(any_of(c("outcome", metrics))) %>%
    pivot_longer(
      cols = -outcome,
      names_to = "metric",
      values_to = "value"
    )

  if (!is.null(outcome_order)) {
    dat$outcome <- factor(dat$outcome, levels = outcome_order)
  }

  ggplot(dat, aes(x = outcome, y = value)) +
    geom_col() +
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw() +
    labs(
      x = NULL,
      y = "Summary value",
      title = "Robustness summary"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_top_influential_subjects <- function(
  loo_all,
  metric = "dfbetas_term",
  top_n = 10,
  use_absolute = TRUE,
  outcome_order = NULL
) {
  dat <- loo_all %>%
    mutate(
      score = if (use_absolute) abs(.data[[metric]]) else .data[[metric]]
    ) %>%
    group_by(outcome) %>%
    slice_max(order_by = score, n = top_n, with_ties = FALSE) %>%
    ungroup()

  if (!is.null(outcome_order)) {
    dat$outcome <- factor(dat$outcome, levels = outcome_order)
  }

  ggplot(dat, aes(x = reorder(subject, score), y = score)) +
    geom_point(size = 2.5) +
    coord_flip() +
    facet_wrap(~ outcome, scales = "free_y") +
    theme_bw() +
    labs(
      x = "Subject",
      y = if (use_absolute) paste0("|", metric, "|") else metric,
      title = paste("Top", top_n, "most influential subjects")
    )
}

plot_metric_rank_across_outcomes <- function(
  metrics_all,
  metric = "R2_marginal"
) {
  dat <- metrics_all %>%
    select(outcome, all_of(metric)) %>%
    arrange(.data[[metric]])

  ggplot(dat, aes(x = reorder(outcome, .data[[metric]]), y = .data[[metric]])) +
    geom_point(size = 3) +
    coord_flip() +
    theme_bw() +
    labs(
      x = NULL,
      y = metric,
      title = paste("Outcome ranking by", metric)
    )
}

plot_influential_counts <- function(
  loo_all,
  outcome_order = NULL
) {
  dat <- loo_all %>%
    group_by(outcome) %>%
    summarise(
      n_influential = sum(influential_dfbetas, na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(outcome_order)) {
    dat$outcome <- factor(dat$outcome, levels = outcome_order)
  }

  ggplot(dat, aes(x = outcome, y = n_influential)) +
    geom_col() +
    theme_bw() +
    labs(
      x = NULL,
      y = "Number of influential subjects",
      title = "Influential subject counts"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

make_all_plots <- function(
  metrics_all,
  loo_all,
  full_metrics = c("R2_marginal", "R2_conditional", "AIC", "RMSE"),
  loo_metrics = c("dfbetas_term", "delta_AIC", "delta_R2_marginal"),
  summary_metrics = c("median_abs_dfbetas", "max_abs_dfbetas", "sd_delta_AIC", "sd_delta_R2_marginal"),
  rank_metric = "R2_marginal",
  top_metric = "dfbetas_term",
  top_n = 10
) {
  outcome_order <- metrics_all %>%
    arrange(.data[[rank_metric]]) %>%
    pull(outcome)

  loo_summary <- summarise_loo_robustness(loo_all)

  plots <- list(
    full_model_metrics = plot_full_model_metrics(
      metrics_all = metrics_all,
      metrics = full_metrics,
      outcome_order = outcome_order
    ),
    loo_distributions = plot_loo_metric_distributions(
      loo_all = loo_all,
      metrics = loo_metrics,
      outcome_order = outcome_order
    ),
    top_influential_subjects = plot_top_influential_subjects(
      loo_all = loo_all,
      metric = top_metric,
      top_n = top_n,
      use_absolute = TRUE,
      outcome_order = outcome_order
    ),
    robustness_summary = plot_loo_robustness_summary(
      loo_summary = loo_summary,
      metrics = summary_metrics,
      outcome_order = outcome_order
    ),
    outcome_ranking = plot_metric_rank_across_outcomes(
      metrics_all = metrics_all,
      metric = rank_metric
    ),
    influential_counts = plot_influential_counts(
      loo_all = loo_all,
      outcome_order = outcome_order
    ),
    loo_summary_table = loo_summary
  )

  return(plots)
}