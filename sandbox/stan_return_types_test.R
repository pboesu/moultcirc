y = as.numeric(circular::rvonmises(n=20, mu = pi, kappa = 5, control.circular = list(modulo = 'asis', zero = 0))) - pi

hist(y)

simple_model <- cmdstanr::cmdstan_model('inst/cmdstanr/simple_vM.stan')

fit <- simple_model$sample(data = list(y = y, N = length(y)))

sm <- fit$summary()
