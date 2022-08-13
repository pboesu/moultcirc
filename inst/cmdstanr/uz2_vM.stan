//
// This Stan program defines the circular Underhill-Zucchini Type 2 model
// The parameterisation follows the moult R package, which may or may not be identical to the parametrisation given in Eq. 5 of Underhill & Zuchini 1988
//

// The input data is a vector 'y' of length 'N'.
data {
  //responses
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_old+N_moult+N_new,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_old+N_moult+N_new,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_kappa;//number of predictors for start date sigma
  matrix[N_old+N_moult+N_new,N_pred_kappa] X_kappa;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=-1*pi(),upper=pi()> alpha_mu;//mean start date // model seems to fail when this is bounded on 0,2*pi
  real<lower=0,upper=2*pi()> alpha_tau;//moult duration
  real<lower=0> alpha_kappa;//von Mises concentration parameter of moult start date
  vector[N_pred_mu-1] beta_mu;//regression coefficients for start date
  vector[N_pred_tau-1] beta_tau;//regression coefficients for duration
  vector[N_pred_kappa-1] beta_kappa;//regression coefficients for sigma start



}

transformed parameters{
  vector[N_old+N_moult+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_new] kappa;//duration lin pred
   mu = X_mu * append_row(alpha_mu,beta_mu);
  //  print(mu);
    tau = X_tau * append_row(alpha_tau,beta_tau);
  //  print(tau);
    kappa = exp(X_kappa * append_row(alpha_kappa,beta_kappa));//use log link for dispersion lin pred
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N_old] P;
  vector[N_moult] q;
  vector[N_new] R;

for (i in 1:N_old) P[i] = 1 - von_mises_cdf(old_dates[i] | mu[i], kappa[i]);
for (i in 1:N_moult) {
  if ((moult_dates[i] - moult_indices[i]*tau[i+N_old]) > -1*pi()) {
    q[i] = log(tau[i+N_old]) + von_mises_lpdf((moult_dates[i] - moult_indices[i]*tau[i+N_old]) | mu[i+N_old], kappa[i+N_old]);//
  } else {
    q[i] = log(tau[i+N_old]) + von_mises_lpdf(((moult_dates[i] - moult_indices[i]*tau[i+N_old])+ 2*pi()) | mu[i+N_old], kappa[i+N_old]);//
  }
}
for (i in 1:N_new) {
  if ((new_dates[i] - tau[i+N_old+N_moult]) > -1*pi()) {
    R[i] = von_mises_cdf((new_dates[i] - tau[i+N_old+N_moult])| mu[i+N_old+N_moult],kappa[i+N_old+N_moult]);//this needs tweaking to ensure the difference is always in [-pi,pi]
  } else {
    R[i] = von_mises_cdf(((new_dates[i] - tau[i+N_old+N_moult]) + 2*pi())| mu[i+N_old+N_moult],kappa[i+N_old+N_moult]);//this needs tweaking to ensure the difference is always in [-pi,pi]
  }
}

target += sum(log(P))+sum(q)+sum(log(R));
//priors
alpha_mu ~ von_mises(0,0.5);
alpha_tau ~ normal(0,pi());//should be truncated? or are the parameter constraints sufficient?
alpha_kappa ~ normal(0,10);
}

generated quantities{
real end_date;
real sigma_approx;
real sigma_exact;
real mu_days;
real tau_days;
 mu_days = (alpha_mu + pi())*365/(2*pi());
 tau_days = alpha_tau*365/(2*pi());
 end_date = mu_days + tau_days;
 sigma_approx = sqrt(1/alpha_kappa)*365/(2*pi());
 sigma_exact = sqrt(1 - modified_bessel_first_kind(1,alpha_kappa)/modified_bessel_first_kind(0,alpha_kappa))*365/(2*pi());//the circular standard deviation -- double check this, as this uses sqrt(circular variance) which is possibly incorrect
}

