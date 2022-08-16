#' Bayesian inference for the Circular Type 2 moult model
#'
#' @export
#' @param date_column the name the column in `data` containing sampling dates, encoded as days since an arbitrary reference date, i.e. a numeric vector
#' @param moult_index_column the name the column in `data` containing moult indices, i.e. a numeric vector of (linearized) moult scores (0 = old plumage,1 = new plumage).
#' @param start_formula model formula for start date
#' @param duration_formula model formula for duration
#' @param kappa_formula model formula for start date sigma
#' @param lump_non_moult logical; should pre- and post-moult observations be treated as indistinguishable? if TRUE, the type 2L model will be fitted.
#' @param data Input data frame must contain a numeric column "date" and a column "moult_index" which is a numeric vector of moult scores ranging from 0 (old plumage) to  1 (new plumage).
#' @param init Specification of initial values for all or some parameters. Can be the string "auto" for an automatic guess based on the data, or any of the permitted rstan options: the digit 0, the strings "0" or "random", or a function. See the detailed documentation for the init argument in ?rstan::stan.
#' @param log_lik boolean retain pointwise log-likelihood in output? This enables model assessment and selection via the loo package. Defaults to true, can lead to very large output arrays if sample size is large.
#' @param ... Arguments passed to `cmdstanr` sampling method (e.g. chains).
#' @importFrom stats model.matrix
#' @return An object of class `moultmcmc`
#'
#TODO: implement an input data class which ensures column names and correct encoding for categorical variables
uz2_circ <- function(moult_index_column,
                     date_column,
                     start_formula = ~1,
                     duration_formula = ~1,
                     kappa_formula = ~1,
                     lump_non_moult = FALSE,
                     data,
                     init = "auto",
                     log_lik = TRUE,
                     ...) {
  stopifnot(all(data[[moult_index_column]] >= 0 & data[[moult_index_column]] <= 1))
  stopifnot(is.numeric(data[[date_column]]))
  stopifnot(is.data.frame(data))
  stopifnot(all(data[[date_column]] >= 0 & data[[date_column]] <= 365))
  #transform dates to radians [-pi,pi]
  data$circday__ = data[[date_column]]*(2*pi / 365 ) - pi
  #order data by moult category
  data <- data[order(data[[moult_index_column]]),]
  #setup model matrices
  X_mu <- model.matrix(start_formula, data)
  X_tau <- model.matrix(duration_formula, data)
  X_kappa <- model.matrix(kappa_formula, data)
  #prepare data structure for stan
  standata <- list(old_dates = data[["circday__"]][data[[moult_index_column]]==0],
                   N_old = length(data[["circday__"]][data[[moult_index_column]]==0]),
                   moult_dates  = data[["circday__"]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)],
                   moult_indices = data[[moult_index_column]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)],
                   N_moult = length(data[["circday__"]][(data[[moult_index_column]] > 0 & data[[moult_index_column]] < 1)]),
                   new_dates = data[["circday__"]][data[[moult_index_column]]==1],
                   N_new = length(data[["circday__"]][data[[moult_index_column]]==1]),
                   X_mu = X_mu,
                   N_pred_mu = ncol(X_mu),
                   X_tau = X_tau,
                   N_pred_tau = ncol(X_tau),
                   X_kappa = X_kappa,
                   N_pred_kappa = ncol(X_kappa),
                   lumped = as.numeric(lump_non_moult))
  #include pointwise log_lik matrix  in output?
  if(log_lik){
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept', 'log_lik')
  } else {
    outpars <- c('beta_mu','beta_tau','beta_sigma', 'sigma_intercept')
  }
  #guess initial values
  if(inherits(init, 'character')){
    if (init == "auto"){
      moult_dates_circ <- circular::circular(standata$moult_dates)
    mu_start = as.numeric(circular::mean.circular(moult_dates_circ))#use circular mean of moult obs
    tau_start = 4*sqrt(1/circular::mle.vonmises(moult_dates_circ)$kappa)#use 4sd of gaussian approximation of dispersion of moult dates
    kappa_start = min(0.5, circular::mle.vonmises(moult_dates_circ)$kappa)# vM MLE of kappa seems to be too tight and then leads to slow sampling, try large dispersion first
    initfunc <- function(chain_id = 1) {
      # cat("chain_id =", chain_id, "\n")
      list(alpha_mu = mu_start, #initialize intercept term from data, set inits for all other effects to 0
           alpha_tau = tau_start,
           alpha_kappa = log(kappa_start),
           beta_mu = as.array(c(rep(0, standata$N_pred_mu - 1))), #initialize intercept term from data, set inits for all other effects to 0
           beta_tau = as.array(c(rep(0, standata$N_pred_tau - 1))),
           beta_kappa = as.array(c(rep(0, standata$N_pred_kappa - 1))))#NB this is on log link scale
    }
    uz2_vM <- cmdstanr::cmdstan_model(system.file("cmdstanr/uz2_vM.stan",package="moultcirc"), stanc_options = list("O1"))
    message(paste('sampling with automatically generated initial values: ', paste(names(unlist(initfunc())), unlist(initfunc()), sep = '=', collapse = ' ')))
    fit <- uz2_vM$sample(data = standata, init = initfunc, ...)
    #out <- rstan::sampling(stanmodels$uz2_linpred, data = standata, init = initfunc, pars = outpars, ...)
  } else { stop('unknown init code')
  }} else {
    uz2_vM <- cmdstanr::cmdstan_model(system.file("cmdstanr/uz2_vM.stan",package="moultcirc"), stanc_options = list("O1"))
    fit <- uz2_vM$sample(data = standata, init = init, ...)
  }
  stanfit <- rstan::read_stan_csv(fit$output_files())
  #rename regression coefficients for output
  # names(out)[grep('beta_mu', names(out))] <- paste('mean',colnames(X_mu), sep = '_')
  # names(out)[grep('beta_tau', names(out))] <- paste('duration',colnames(X_tau), sep = '_')
  # names(out)[grep('beta_sigma', names(out))] <- paste('log_sd',colnames(X_sigma), sep = '_')
  # names(out)[grep('sigma_intercept', names(out))] <- 'sd_(Intercept)'
  out_struc <- list()
  out_struc$stanfit <- stanfit
  out_struc$terms$date_column <- date_column
  out_struc$terms$moult_index_column <- moult_index_column
  out_struc$terms$moult_cat_column <- NA
  out_struc$stan_profile <- fit$profiles()
  class(out_struc) <- 'moultmcmc'
  return(out_struc)
}
