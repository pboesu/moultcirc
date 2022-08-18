data("sim_data_small")

test_that("model runs through for moult contiguous in [-pi,pi]", {

  uz2_circ_ufit <- uz2_circ("moult_score", "yday", data = sim_data_small[1:50,], chains = 2, parallel_chains = 2, refresh = 500)
  expect_s3_class(uz2_circ_ufit, "moultmcmc")
  #get quantiles
  q_days =quantile(rstan::extract(uz2_circ_ufit$stanfit,pars = "mu_days")$mu_days, p = c(0.025, 0.975))
  q_duration = quantile(rstan::extract(uz2_circ_ufit$stanfit,pars = "tau_days")$tau_days, p = c(0.025, 0.975))
  testthat::expect_equal(TRUE, dplyr::between(157, q_days[1],q_days[2]))
  testthat::expect_equal(TRUE, dplyr::between(102, q_duration[1],q_duration[2]))
  testthat::expect_equal(TRUE, all(rstan::summary(uz2_circ_ufit$stanfit)$summary[,'Rhat']<1.1))

  })

test_that("model runs through for moult data that are bimodal in [-pi,pi] ", {

  uz2_circ_usfit <- uz2_circ("moult_score", "yday_shifted", data = sim_data_small[1:50,], chains = 2, parallel_chains = 2, refresh = 500)
  expect_s3_class(uz2_circ_usfit, "moultmcmc")
  #get quantiles
  q_days =quantile(rstan::extract(uz2_circ_usfit$stanfit,pars = "mu_days")$mu_days, p = c(0.025, 0.975))
  q_duration = quantile(rstan::extract(uz2_circ_usfit$stanfit,pars = "tau_days")$tau_days, p = c(0.025, 0.975))
  testthat::expect_equal(TRUE, dplyr::between(157+150, q_days[1],q_days[2]))
  testthat::expect_equal(TRUE, dplyr::between(102, q_duration[1],q_duration[2]))
  testthat::expect_equal(TRUE, all(rstan::summary(uz2_circ_usfit$stanfit)$summary[,'Rhat']<1.1))
})

test_that("lumped model runs through for moult in [-pi,pi]", {

  uz2_circ_lfit <- uz2_circ("moult_score", "yday", lump_non_moult = TRUE, data = sim_data_small[1:50,], refresh = 100, iter_warmup = 100, iter_sampling = 100, chains = 2, parallel_chains = 2)#very slow 35mins for 300iter
  expect_s3_class(uz2_circ_lfit, "moultmcmc")
  #get quantiles
  q_days =quantile(rstan::extract(uz2_circ_lfit$stanfit,pars = "mu_days")$mu_days, p = c(0.025, 0.975))
  q_duration = quantile(rstan::extract(uz2_circ_lfit$stanfit,pars = "tau_days")$tau_days, p = c(0.025, 0.975))
  message(paste(q_days, 'days'))
  message(paste(q_duration, 'days'))
  testthat::expect_equal(TRUE, dplyr::between(157, q_days[1],q_days[2]))
  testthat::expect_equal(TRUE, dplyr::between(102, q_duration[1],q_duration[2]))
  testthat::expect_equal(TRUE, all(rstan::summary(uz2_circ_lfit$stanfit)$summary[,'Rhat']<1.1))

})
