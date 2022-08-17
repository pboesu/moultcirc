//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as von mises distributed
// with mean 'mu' and concentration 'kappa'.
//
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector<lower=-1*pi(),upper=pi()>[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=-1*pi(),upper=pi()> mu;
  real<lower=0> kappa;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ von_mises(mu, kappa);
}

generated quantities{
  real P;
  vector[N] PvM;
  real ProdvM;
  real PN;
  real SumPhiN;
  real ProdN;
  vector[N] PhiN;

  //real random;
  //array[2] real random2;

  P = von_mises_cdf(y | mu, kappa);
  for (i in 1:N) PvM[i] = von_mises_cdf(y[i] | mu, kappa);
  ProdvM = prod(PvM);
  //P = von_mises_lpdf(y | mu, kappa);
  //P = von_mises_lupdf(y | mu, kappa);
  //P = von_mises_lcdf(y | mu, kappa);
  //P = von_mises_lccdf(y | mu, kappa);
  PN = normal_cdf(y | mu, sqrt(1/kappa));
  PhiN = Phi((y-mu)/sqrt(1/kappa));
  SumPhiN = sum(PhiN);
  ProdN = prod(PhiN);
  //random = von_mises_rng(mu, kappa);
  //random2 = von_mises_rng(y[1:2], kappa);



}

