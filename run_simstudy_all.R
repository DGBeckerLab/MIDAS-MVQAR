N <- 1000
for (param_set in 1:4) {
  for (theta_shrink in c(0, 2e-3, 1e-2)) {
    source("run_simstudy.R")
  }
}
