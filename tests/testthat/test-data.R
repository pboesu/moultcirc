data("sim_data_small")

test_that("model runs through for moult in [-pi,pi]", {

  uz2_circ_fit <- uz2_circ("moult_score", "yday", data = sim_data_small, chains = 1, refresh = 500)
  expect_s3_class(uz2_circ_fit, "moultmcmc")
  #get quantiles
  q_days =quantile(rstan::extract(uz2_circ_fit$stanfit,pars = "mu_days")$mu_days, p = c(0.025, 0.975))
  q_duration = quantile(rstan::extract(uz2_circ_fit$stanfit,pars = "tau_days")$tau_days, p = c(0.025, 0.975))
  testthat::expect_equal(TRUE, dplyr::between(157, q_days[1],q_days[2]))
  testthat::expect_equal(TRUE, dplyr::between(102, q_duration[1],q_duration[2]))
  })
