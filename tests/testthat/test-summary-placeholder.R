source(testthat::test_path("..", "..", "functions.R"))

test_that("summary keeps tau positions when a CI table is NULL", {
  tab <- data.frame(
    param = paste0("p", 1:42),
    estimate = seq_len(42) / 10,
    se = rep(0.1, 42),
    ci_lower = seq_len(42) / 10 - 0.2,
    ci_upper = seq_len(42) / 10 + 0.2
  )
  fit <- list(
    "tau_0.25" = list(ci_table = NULL, loss = 1, convergence = 0),
    "tau_0.50" = list(ci_table = tab, loss = 2, convergence = 0)
  )
  class(fit) <- "mvqarmidas_fit"
  attr(fit, "p") <- 1
  attr(fit, "taus") <- c(0.25, 0.5)

  s <- summary(fit)
  expect_length(s, 2)
  expect_null(s[[1]])
  expect_equal(s[[2]]$tau, 0.5)
  expect_equal(s[[2]]$loss, 2)
})

test_that("MIDAS block columns are ordered lag0 = most recent month", {
  x <- seq_len(20)
  blocks <- create_midas_blocks(x, n = 6, m = 3, K = 4, prefix = "g")
  # Row 6 ends at month 15, so lag0 should be 15 and lag4 should be 11
  expect_equal(unname(blocks[6, ]), c(15, 14, 13, 12, 11))
  expect_equal(colnames(blocks), paste0("g_lag", 0:4))
})
