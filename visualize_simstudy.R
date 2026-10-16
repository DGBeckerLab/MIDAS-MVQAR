library(dplyr)
library(tidyr)
library(ggplot2)

N_plot <- 250
source("functions.R")

files <- list.files(
  path = "results",
  pattern = sprintf(
    "^sim-m3-K4-p1-tau0\\.5-param_set[0-9]+-lambda[0-9.e-]+-N%d\\.csv$",
    N_plot
  ),
  full.names = TRUE
)
res <- bind_rows(lapply(files, read.csv))
count(res, param_set, theta_shrink, name = "rows")

truth <- bind_rows(lapply(1:4, function(k) {
  tp <- create_true_par_vector(1, param_set = k)
  data.frame(param_set = k, param = names(tp), true = unname(tp))
}))

par_names <- grep("^eq", names(res), value = TRUE)
cov_names <- paste0("cov_", par_names)
ok <- filter(res, status == "ok")

est <- ok |>
  select(param_set, theta_shrink, sim, all_of(par_names)) |>
  pivot_longer(all_of(par_names), names_to = "param", values_to = "est")

cvr <- ok |>
  select(param_set, theta_shrink, sim, all_of(cov_names)) |>
  pivot_longer(all_of(cov_names), names_to = "param", values_to = "covered") |>
  mutate(param = sub("^cov_", "", param))

d <- est |>
  left_join(cvr, by = c("param_set", "theta_shrink", "sim", "param")) |>
  left_join(truth, by = c("param_set", "param")) |>
  mutate(err = est - true)

summ <- d |>
  group_by(param_set, theta_shrink, param, true) |>
  summarise(
    n_ok = n(),
    bias = mean(err),
    se_bias = sd(est) / sqrt(n_ok),
    rmse = sqrt(mean(err^2)),
    coverage = mean(covered, na.rm = TRUE),
    se_cov = sqrt(0.95 * 0.05 / n_ok),
    .groups = "drop"
  ) |>
  mutate(
    eq = sub("_.*", "", param),
    # Parameter subsets are plotted separately
    group = case_when(
      grepl("alpha", param) ~ "alpha",
      grepl("midas_scale", param) ~ "midas_scale",
      grepl("theta", param) ~ "theta",
      grepl("phi", param) ~ "phi"
    ),
    param = factor(param, levels = par_names)
  )

status_tab <- res |>
  mutate(status_type = ifelse(startsWith(status, "error"), "error", status)) |>
  count(param_set, theta_shrink, status_type)

plot_metric <- function(dat, value, xlab, se = NULL, ref = NULL) {
  p <- ggplot(dat, aes(y = param, x = .data[[value]], color = factor(theta_shrink))) +
    facet_grid(eq ~ param_set, labeller = label_both, scales = "free_y") +
    labs(x = xlab, y = NULL, color = "lambda")
  if (!is.null(ref)) p <- p + geom_vline(xintercept = ref, linetype = "dashed")
  if (is.null(se)) {
    p + geom_point(position = position_dodge(width = 0.6))
  } else {
    p + geom_pointrange(
      aes(xmin = .data[[value]] - 2 * .data[[se]], xmax = .data[[value]] + 2 * .data[[se]]),
      position = position_dodge(width = 0.6)
    )
  }
}

dir.create("figures", showWarnings = FALSE)

for (g in unique(summ$group)) {
  dg <- filter(summ, group == g)
  ggsave(sprintf("figures/bias_%s_N%d.png", g, N_plot),
    plot_metric(dg, "bias", "Bias (±2 MC SE)", se = "se_bias", ref = 0),
    width = 10, height = 9, dpi = 150, bg = "white")
  ggsave(sprintf("figures/rmse_%s_N%d.png", g, N_plot),
    plot_metric(dg, "rmse", "RMSE"),
    width = 10, height = 9, dpi = 150, bg = "white")
  ggsave(sprintf("figures/coverage_%s_N%d.png", g, N_plot),
    plot_metric(dg, "coverage", "95% CI coverage (±2 SE)", se = "se_cov", ref = 0.95),
    width = 10, height = 9, dpi = 150, bg = "white")
}

ggsave(sprintf("figures/status_N%d.png", N_plot),
  ggplot(status_tab, aes(x = factor(theta_shrink), y = n, fill = status_type)) +
    geom_col(position = "fill") +
    facet_wrap(~param_set, labeller = label_both) +
    labs(x = "lambda", y = "share of replications", fill = "status"),
  width = 8, height = 5, dpi = 150, bg = "white")
