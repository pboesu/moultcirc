library(circular)
dvonmises()
plot(pvonmises(seq(0,2*pi, by =0.1), mu = pi, kappa = 5))
plot(1-pvonmises(seq(0,2*pi, by =0.1), mu = pi, kappa = 5))
points(pvonmises(seq(0,2*pi, by =0.1)-pi/2, mu = pi, kappa = 5), pch = 16)

plot(1 - pvonmises(seq(0,2*pi, by =0.01), mu = pi, kappa = 5) + (pvonmises(seq(0,2*pi, by =0.01)-pi/2, mu = pi, kappa = 5)))
lines(pmax(1 - pvonmises(seq(0,2*pi, by =0.01), mu = pi, kappa = 5) , (pvonmises(seq(0,2*pi, by =0.01)-pi/2, mu = pi, kappa = 5))))
lines(log(exp(1 - pvonmises(seq(0,2*pi, by =0.01), mu = pi, kappa = 5))+ exp(pvonmises(seq(0,2*pi, by =0.01)-pi/2, mu = pi, kappa = 5))))
#try to this also with the stan functions using exposed stan functions


days = 1:365
mu = 150
tau = 300
sigma = 7

P = 1 - pnorm(days, mu, sigma)
Q = pnorm(days, mu, sigma) - pnorm(days - tau, mu, sigma)
R = pnorm(days - tau, mu, sigma)
plot(days, P)
points(days, Q, col = 'red')
points(days, R, col = 'blue')

notQ = (1 - pnorm(days, mu, sigma)) + pnorm(days - tau, mu, sigma)
lines(days, notQ, col = 'green')


angles = seq(-pi,pi,by = 0.05)
mua = mu/365*2*pi-pi
taua = tau/365*2*pi
kappa = 1/(sigma/365*2*pi)^2

Pa = 1 - pvonmises(angles, mua, kappa)
Qa = pvonmises(angles, mua, kappa) - pvonmises(angles - taua, mua, kappa)
Ra = pvonmises(angles - taua, mua, kappa)
plot(angles, Pa, ylim = c(-1,2))
points(angles, Qa, col = 'red')
points(angles, Ra, col = 'blue')
notQa = 1 - Qa
lines(angles, notQa)
notQa2 = pmax(Pa, Ra)
notQa3 = pmin(1, Pa+Ra)
lines(angles, notQa2, col= 'darkgreen', lty = 2)
lines(angles, notQa3, col= 'darkred', lty = 2, lwd = 3)




library(Rcpp)
Rcpp::cppFunction(
  'double besselJ(double v, double z) {
  return stan::math::modified_bessel_first_kind(v, z);
}', depends = "StanHeaders", includes = "#include <stan/math.hpp>")
