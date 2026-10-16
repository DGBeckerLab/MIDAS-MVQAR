suppressPackageStartupMessages({
  library(xts)
  library(zoo)
  library(quantreg)
  library(midasr)
  library(dplyr)
  library(stringr)
})

source("functions.R")

# param_set, theta_shrink, and N must be set before sourcing (run_simstudy_all.R does this).
# Other settings are fixed here so leftover globals cannot change a run silently.
stopifnot(exists("param_set"), exists("theta_shrink"), exists("N"))
stopifnot(length(param_set) == 1, length(theta_shrink) == 1)
taus <- 0.5
theta_bound <- 5
start_method <- "neutral"
kappa <- 0.01
n_restarts <- 2
ci_method <- "sandwich"
force_rerun <- FALSE

m <- 3
K <- 4
p <- 1
seed <- 6173

dir.create("results", showWarnings = FALSE)

sim_file_name <- paste0(
  "results/sim-m", m, "-K", K, "-p", p, "-tau", taus,
  "-param_set", param_set, "-lambda", sprintf("%g", theta_shrink),
  "-start", start_method, "-tb", theta_bound, "-kappa", kappa,
  "-seed", seed, "-restarts", n_restarts, "-ci", ci_method,
  "-N", N, ".csv"
)
message("Will save to file: ", sim_file_name)

true_par <- create_true_par_vector(p, param_set = param_set)
par_names <- names(true_par)

summarise_cell <- function(res_df, true_par) {
  ok <- res_df$status == "ok"
  est <- as.matrix(res_df[ok, par_names, drop = FALSE])
  cov <- as.matrix(res_df[ok, paste0("cov_", par_names), drop = FALSE])
  n_ok <- nrow(est)
  bias <- colMeans(est) - true_par
  data.frame(
    param = par_names,
    true = true_par,
    mean_est = colMeans(est),
    bias = bias,
    mc_se_bias = apply(est, 2, sd) / sqrt(n_ok),
    rmse = sqrt(colMeans((est - matrix(true_par, n_ok, length(true_par), byrow = TRUE))^2)),
    coverage95 = colMeans(cov, na.rm = TRUE),
    n_ok = n_ok,
    n_cov_na = colSums(is.na(cov)),
    row.names = NULL
  )
}

if (file.exists(sim_file_name) && !force_rerun) {
  res_df <- read.csv(sim_file_name)
} else {
  set.seed(seed)
  start_par <- create_neutral_par_vector(p)

  res_est <- matrix(NA_real_, nrow = N, ncol = length(true_par), dimnames = list(NULL, par_names))
  res_cov <- matrix(NA_real_, nrow = N, ncol = length(true_par), dimnames = list(NULL, paste0("cov_", par_names)))
  status <- character(N)

  t0 <- Sys.time()
  for (i in seq_len(N)) {
    sim_data <- simulate(m = m, K = K, p = p, taus = taus, param_set = param_set)

    fit <- tryCatch(
      estimate_joint_qarmidas(
        taus = taus,
        data = sim_data,
        K = K,
        p = p,
        kappa = kappa,
        theta_shrink = theta_shrink,
        start_par = if (start_method == "ols") create_ols_start_par(sim_data, K, p) else start_par,
        theta_bound = theta_bound,
        optimizer = "L-BFGS-B",
        use_bounds = FALSE,
        n_restarts = n_restarts,
        jitter_sd = 0.01,
        verbose = FALSE,
        dry_run = FALSE,
        use_analytic_grad = TRUE,
        ci_method = ci_method
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      status[i] <- paste0("error: ", conditionMessage(fit))
      next
    }

    tau_fit <- fit[[1]]
    est <- tau_fit$opt_result$par
    res_est[i, ] <- est
    ci <- tau_fit$ci_table
    res_cov[i, ] <- if (is.null(ci)) NA_real_ else as.numeric(ci$ci_lower <= true_par & true_par <= ci$ci_upper)
    status[i] <- if (tau_fit$convergence == 0) "ok" else "maxit"

    if (i %% 25 == 0) {
      cat("Completed", i, "of", N, "after",
          round(difftime(Sys.time(), t0, units = "mins"), 1), "mins\n")
    }
  }

  res_df <- data.frame(
    sim = seq_len(N),
    status = status,
    res_est,
    res_cov,
    row.names = NULL
  )
  res_df <- cbind(
    param_set = param_set,
    theta_shrink = theta_shrink,
    res_df
  )
  write.csv(res_df, file = sim_file_name, row.names = FALSE)
  cat("Completed in", round(difftime(Sys.time(), t0, units = "mins"), 2), "mins\n")
}

cat("\nStatus counts (param_set =", param_set, ", theta_shrink =", theta_shrink, "):\n")
print(table(res_df$status))

summary_tab <- summarise_cell(res_df, true_par)
summary_tab
