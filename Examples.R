# library(devtools)
# install_github("https://github.com/jorgeyslas/phfrailty.git")

library(phfrailty)
library(expm)
library(xtable)
library(cubature)
library(copula)
source("Gamma_model.R")
source("bms_ph_fitting.R")
source("ph_premiums.R")


# Simulated data - Independence
M <- 2500
set.seed(1)
ind_dat <- sim_ind_gamma(M, 10 / 3, 10 / 3, 10 / 3, 10 / 3)

theta1 <- ind_dat$theta1
theta2 <- ind_dat$theta2
n <- ind_dat$n
w <- ind_dat$w
y <- ind_dat$y

### Frequency modeling
p <- 3
set.seed(1)
mph_ini <- mp(phasetype(structure = "gcoxian", dimension = p))
mph_fit <- fit(mph_ini, n, stepsEM = 100, stepsPH = 50, every = 10)

# Loglikelihood using true generating model
sum(log(gamma_mix_densf(n, 10 / 3, 10 / 3)))


pdf("freq_hist.pdf", width = 6, height = 6)
plot(table(n) / M,
  lwd = 2, ylim = c(0, 0.43),
  ylab = "", xlab = "", main = "Histogram vs fitted density \n (Frequency)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5, col = "darkgray"
)
points(0:max(n), gamma_mix_densf(0:max(n), 10 / 3, 10 / 3), col = "#3498DB", lw = 2, pch = 3, cex = 1.2)
points(0:max(n), dens(mph_fit, 0:max(n)), col = "#E74C3C", lw = 2, cex = 1.2)
legend("topright", inset = 0.05, bty = "n", c("Original density", "Fitted density"), horiz = FALSE, pch = c(3, 1), col = c("#3498DB", "#E74C3C"), cex = 1.2)
dev.off()

under_ph <- phasetype(alpha = mph_fit@pars$alpha, S = mph_fit@pars$S)
sq <- seq(0, 5, 0.01)
pdf("freq_under.pdf", width = 6, height = 6)
plot(sq, dgamma(sq, 10 / 3, 10 / 3),
  type = "l", lwd = 2,
  ylab = "", xlab = expression(theta[1]), main = "Original vs fitted mixing densities \n (Frequency)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5
)
lines(sq, dens(under_ph, sq), col = "#E74C3C", lwd = 2, lty = 5)
legend("topright", inset = 0.05, bty = "n", c("Original mixing density", "Fitted mixing density"), lwd = 2, lty = c(1, 5), col = c(1, 2), cex = 1.2)
dev.off()

# Moments
mean(n)
moment(under_ph, 1)

# Fitted parameters
round(mph_fit@pars$alpha, digits = 4)
round(mph_fit@pars$S, digits = 4)


### Severity modeling
wp <- w[n > 0]
np <- n[n > 0]

p <- 3
set.seed(1)
ph_ini <- phasetype(structure = "gcoxian", dimension = p)

ph_fit <- fit_ph(ph_ini, np, wp, stepsEM = 100, stepsPH = 50, every = 1)
fph_fit <- frailty(ph_fit, bhaz = "exponential", bhaz_pars = 1)

# Loglikelihood using true generating model
sum(log(gamma_mix_denssl(np, wp, 10 / 3, 10 / 3)))

# Plot
sq <- seq(0, 8, by = 0.01)
pdf("sev_cdf.pdf", width = 6, height = 6)
Fhat <- ecdf(y)
plot(Fhat,
  verticals = TRUE, do.points = FALSE,
  xlab = "", ylab = "",
  main = "Empirical CDF vs fitted CDF \n (Severity)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5,
  xlim = range(sq), ylim = c(0, 1), lwd = 2
)
lines(sq, gamma_mix_cdf(sq, 10 / 3, 10 / 3), col = "#3498DB", lwd = 2)
lines(sq, cdf(fph_fit, sq), col = "#E74C3C", lwd = 2, lty = 5)
legend("bottomright",
  inset = 0.05, bty = "n",
  legend = c("Empirical CDF", "Original CDF", "Fitted CDF"),
  lty = c(1, 1, 5), lwd = 2,
  col = c("black", "#3498DB", "#E74C3C"),
  cex = 1.2
)
dev.off()

pdf("sev_under.pdf", width = 6, height = 6)
sq <- seq(0, 5, 0.01)
plot(sq, dgamma(sq, 10 / 3, 10 / 3),
  type = "l", lwd = 2,
  ylab = "", xlab = expression(theta[1]), main = "Original vs fitted mixing densities \n (Severity)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5
)
lines(sq, dens(ph_fit, sq), col = "#E74C3C", lwd = 2, lty = 5)
legend("topright", inset = 0.05, bty = "n", c("Original mixing density", "Fitted mixing density"), lwd = 2, lty = c(1, 5), col = c(1, 2), cex = 1.2)
dev.off()

# Moments
sum(w) / M
mean(y)

inverse_ph_moment(ph_fit)
round(ph_fit@pars$alpha %*% (-ph_fit@pars$S) %*% rep(1, p), 6)
integrate(function(x) dgamma(x, 10 / 3, 10 / 3) / x, 0, Inf)$value

# Fitted Parameters
round(ph_fit@pars$alpha, digits = 4)
round(ph_fit@pars$S, digits = 4)


### Premiums

# Gamma model
years <- 1:5
nclaims <- 1:5
wo <- 3.5
premiums <- prem_gamma_mat(nclaims, years, wo, 10 / 3, 10 / 3, 10 / 3, 10 / 3)

round(premiums, 4)
rownames(premiums) <- 0:5

print(
  xtable(premiums,
    digits = 4,
    caption = "Optimal net premiums under the true Gamma model."
  ),
  include.rownames = TRUE
)

# Independent PH model
premiums_ph <- prem_ph_mat(nclaims, years, wo, mph_fit, ph_fit)
round(premiums_ph, 4)
rownames(premiums_ph) <- 0:5
print(
  xtable(premiums_ph,
    digits = 4,
    caption = "Optimal net premiums under the ind PH model."
  ),
  include.rownames = TRUE
)


### Joint frequency-severity modeling
p1 <- 3
p2 <- 3
set.seed(123)
bivph_ini <- bivphasetype(dimensions = c(p1, p2))

bmph_fit <- fit_bph(
  bivph_ini,
  n,
  w,
  stepsEM = 100,
  stepsPH = 100,
  nq1 = 200,
  nq2 = 200,
  upper1 = 8,
  upper2 = 8,
  keep_mass = 0.9995,
  keep_ratio = 0.9995,
  every = 1
)
fit_biv <- bmph_fit$bph_fit

# Loglikelihood using true underlying model
ind_gamma(n, w, 10 / 3, 10 / 3, 10 / 3, 10 / 3)


ph1 <- marginal(fit_biv, 1)
fbph_fit <- mp(ph1)

pdf("freq_hist_biv.pdf", width = 6, height = 6)
plot(table(n) / M,
  lwd = 2, ylim = c(0, 0.43),
  ylab = "", xlab = "", main = "Histogram vs fitted density \n (Frequency - Joint modeling)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5, col = "darkgray"
)
points(0:max(n), gamma_mix_densf(0:max(n), 10 / 3, 10 / 3), col = "#3498DB", lw = 2, pch = 3, cex = 1.2)
points(0:max(n), dens(mph_fit, 0:max(n)), col = "#E74C3C", lw = 2, cex = 1.2)
legend("topright", inset = 0.05, bty = "n", c("Original density", "Fitted density"), horiz = FALSE, pch = c(3, 1), col = c("#3498DB", "#E74C3C"), cex = 1.2)
dev.off()

ph2 <- marginal(fit_biv, 2)
sbph_fit <- frailty(ph2, bhaz = "exponential", bhaz_pars = 1)

sq <- seq(0, 8, by = 0.01)
pdf("sev_cdf_biv.pdf", width = 6, height = 6)
Fhat <- ecdf(y)
plot(Fhat,
  verticals = TRUE, do.points = FALSE,
  xlab = "", ylab = "",
  main = "Empirical CDF vs fitted CDF \n (Severity - Joint modeling)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5,
  xlim = range(sq), ylim = c(0, 1), lwd = 2
)
lines(sq, gamma_mix_cdf(sq, 10 / 3, 10 / 3), col = "#3498DB", lwd = 2)
lines(sq, cdf(sbph_fit, sq), col = "#E74C3C", lwd = 2, lty = 5)
legend("bottomright",
  inset = 0.05, bty = "n",
  legend = c("Empirical CDF", "Original CDF", "Fitted CDF"),
  lty = c(1, 1, 5), lwd = 2,
  col = c("black", "#3498DB", "#E74C3C"),
  cex = 1.2
)
dev.off()


# Moments
ph2@pars$alpha %*% (-ph2@pars$S) %*% rep(1, p2)
biv_ph_moments(fit_biv)
corr(fit_biv)

# Parameters
round(fit_biv@pars$alpha, 4)
round(fit_biv@pars$S11, 4)
round(fit_biv@pars$S12, 4)
round(fit_biv@pars$S22, 4)


# Premiums
premiums_bivph <- prem_bivph_mat(nclaims, years, wo, fit_biv)
round(premiums_bivph, 4)
rownames(premiums_bivph) <- 0:5
print(
  xtable(premiums_bivph,
    digits = 4,
    caption = "Optimal net premiums under the biv PH model."
  ),
  include.rownames = TRUE
)


### Dependent frequency-severity example

## Simulation - rho = 0.5
rho <- 0.5
M <- 2500
set.seed(1)
dep_dat <- sim_dep_gamma(M, rho, 2, 1 / 2, 10 / 3, 10 / 3)

theta1 <- dep_dat$theta1
theta2 <- dep_dat$theta2
n <- dep_dat$n
w <- dep_dat$w
y <- dep_dat$y

round(mean(n == 0), 4)
plot(theta1, theta2)
plot(theta1, 1 / theta2)
cor(theta1, theta2)
cor(theta1, 1 / theta2)

biv_gam_moments(rho, 2, 1 / 2, 10 / 3, 10 / 3)

M_prem <- 100000
set.seed(1)
prem_theta <- sim_latent_dep_gamma(M_prem, rho, 2, 1 / 2, 10 / 3, 10 / 3)

## Independent fitting
# Frequency
p <- 2
set.seed(123)
mph_ini <- mp(phasetype(structure = "gcoxian", dimension = p))
mph_fit <- fit(mph_ini, n, stepsEM = 10, stepsPH = 50, every = 10, truncationpoint = 20)
#-6136.953

# Loglikelihood using true marginal model
sum(log(gamma_mix_densf(n, 2, 1 / 2)))

under_ph <- phasetype(alpha = mph_fit@pars$alpha, S = mph_fit@pars$S)

mean(n)
moment(under_ph, 1)

# Severity
wp <- w[n > 0]
np <- n[n > 0]

p <- 3
set.seed(1)
ph_ini <- phasetype(structure = "gcoxian", dimension = p)

ph_fit <- fit_ph(ph_ini, np, wp, stepsEM = 100, stepsPH = 50, every = 1)
#-10518.57
fph_fit <- frailty(ph_fit, bhaz = "exponential", bhaz_pars = 1)

# Loglikelihood using true marginal model
sum(log(gamma_mix_denssl(np, wp, 10 / 3, 10 / 3)))


sum(w) / M
mean(y)

ph_fit@pars$alpha %*% (-ph_fit@pars$S) %*% rep(1, p)
inverse_ph_moment(ph_fit)

inverse_ph_moment(ph_fit) * moment(under_ph, 1)

## Joint frequency-severity modeling
p1 <- 6
p2 <- 6
set.seed(19)
bivph_ini <- bivphasetype(dimensions = c(p1, p2))

bmph_fit <- fit_bph(
  bivph_ini,
  n,
  w,
  stepsEM = 100,
  stepsPH = 100,
  nq1 = 200,
  nq2 = 200,
  upper1 = 25,
  upper2 = 8,
  keep_mass = 0.9995,
  keep_ratio = 0.9995,
  every = 1
)

fit_biv <- bmph_fit$bph_fit

bmph_fit$loglik_trace[length(bmph_fit$loglik_trace)]

ph1 <- marginal(fit_biv, 1)
fbph_fit <- mp(ph1)

# Frequency plot
pdf("freq_hist_dep_p05.pdf", width = 6, height = 6)
plot(table(n) / M,
  lwd = 2, ylim = c(0, 0.17),
  ylab = "", xlab = "", main = "Histogram vs fitted density \n (Frequency - Rho = 0.5)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5, col = "darkgray"
)
points(0:max(n), gamma_mix_densf(0:max(n), 2, 0.5), col = "#3498DB", lw = 2, pch = 3, cex = 1.2)
points(0:max(n), dens(mph_fit, 0:max(n)), col = "#339966", lw = 2, cex = 1.2, pch = 2)
points(0:max(n), dens(fbph_fit, 0:max(n)), col = "#E74C3C", lw = 2, cex = 1.2)
legend("topright", inset = 0.05, bty = "n", c("Original density", "Independent PH", "Bivariate PH"), horiz = FALSE, pch = c(3, 2, 1), col = c("#3498DB", "#339966", "#E74C3C"), cex = 1.2)
dev.off()

ph2 <- marginal(fit_biv, 2)
sbph_fit <- frailty(ph2, bhaz = "exponential", bhaz_pars = 1)

# Severity plot
sq <- seq(0, 8, by = 0.01)
pdf("sev_cdf_dep_p05.pdf", width = 6, height = 6)
Fhat <- ecdf(y)
plot(Fhat,
  verticals = TRUE, do.points = FALSE,
  xlab = "", ylab = "",
  main = "Empirical CDF vs fitted CDF \n (Severity - Rho = 0.5)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5,
  xlim = range(sq), ylim = c(0, 1), lwd = 2
)
lines(sq, gamma_mix_cdf(sq, 10 / 3, 10 / 3), col = "#3498DB", lwd = 2)
lines(sq, cdf(fph_fit, sq), col = "#339966", lwd = 2, lty = 3)
lines(sq, cdf(sbph_fit, sq), col = "#E74C3C", lwd = 2, lty = 5)
legend("bottomright",
  inset = 0.05, bty = "n",
  legend = c("Empirical CDF", "Original CDF", "Independent PH", "Bivariate PH"),
  lty = c(1, 1, 3, 5), lwd = 2,
  col = c("black", "#3498DB", "#339966", "#E74C3C"),
  cex = 1.2
)
dev.off()


# Moments
ph2@pars$alpha %*% (-ph2@pars$S) %*% rep(1, p2)
biv_ph_moments(fit_biv)
round(corr(fit_biv), 4)

# Premiums
premiums_ind_ph <- prem_ph_mat(nclaims, years, wo, mph_fit, ph_fit)
premiums_biv_ph <- prem_bivph_mat(nclaims, years, wo, fit_biv)
premiums_gamma <- prem_gamma_dep_mat(nclaims, years, wo, prem_theta$theta1, prem_theta$theta2)

idx_updated <- 2:nrow(premiums_gamma)

P_true <- premiums_gamma[idx_updated, ]
P_ind <- premiums_ind_ph[idx_updated, ]
P_biv <- premiums_biv_ph[idx_updated, ]

RMSE_ind <- sqrt(mean((P_ind - P_true)^2))
RMSE_biv <- sqrt(mean((P_biv - P_true)^2))

RRMSE_ind <- sqrt(mean(((P_ind - P_true) / P_true)^2))
RRMSE_biv <- sqrt(mean(((P_biv - P_true) / P_true)^2))

round(RRMSE_ind * 100, 2)
round(RRMSE_biv * 100, 2)


## Simulation - rho = 0.25
rho <- 0.25
M <- 2500
set.seed(1)
dep_dat <- sim_dep_gamma(M, rho, 2, 1 / 2, 10 / 3, 10 / 3)

theta1 <- dep_dat$theta1
theta2 <- dep_dat$theta2
n <- dep_dat$n
w <- dep_dat$w
y <- dep_dat$y

round(mean(n == 0), 4)

biv_gam_moments(rho, 2, 1 / 2, 10 / 3, 10 / 3)

M_prem <- 100000
set.seed(1)
prem_theta <- sim_latent_dep_gamma(M_prem, rho, 2, 1 / 2, 10 / 3, 10 / 3)

## Independent fitting
# Frequency
p <- 2
set.seed(123)
mph_ini <- mp(phasetype(structure = "gcoxian", dimension = p))
mph_fit <- fit(mph_ini, n, stepsEM = 10, stepsPH = 50, every = 10, truncationpoint = 20)
#-6142.31

# Severity
wp <- w[n > 0]
np <- n[n > 0]

p <- 3
set.seed(1)
ph_ini <- phasetype(structure = "gcoxian", dimension = p)

ph_fit <- fit_ph(ph_ini, np, wp, stepsEM = 100, stepsPH = 50, every = 1)
#-11275.94

ph_fit@pars$alpha %*% (-ph_fit@pars$S) %*% rep(1, p)


## Joint frequency-severity modeling
p1 <- 6
p2 <- 6
set.seed(1313)
bivph_ini <- bivphasetype(dimensions = c(p1, p2))

bmph_fit <- fit_bph(
  bivph_ini,
  n,
  w,
  stepsEM = 100,
  stepsPH = 100,
  nq1 = 200,
  nq2 = 200,
  upper1 = 25,
  upper2 = 8,
  keep_mass = 0.9995,
  keep_ratio = 0.9995,
  every = 1
)

fit_biv <- bmph_fit$bph_fit

bmph_fit$loglik_trace[length(bmph_fit$loglik_trace)]

ph1 <- marginal(fit_biv, 1)
fbph_fit <- mp(ph1)

ph2 <- marginal(fit_biv, 2)
sbph_fit <- frailty(ph2, bhaz = "exponential", bhaz_pars = 1)


# Moments
ph2@pars$alpha %*% (-ph2@pars$S) %*% rep(1, p2)
biv_ph_moments(fit_biv)
round(corr(fit_biv), 4)

# Premiums
premiums_ind_ph <- prem_ph_mat(nclaims, years, wo, mph_fit, ph_fit)
premiums_biv_ph <- prem_bivph_mat(nclaims, years, wo, fit_biv)
premiums_gamma <- prem_gamma_dep_mat(nclaims, years, wo, prem_theta$theta1, prem_theta$theta2)

idx_updated <- 2:nrow(premiums_gamma)

P_true <- premiums_gamma[idx_updated, ]
P_ind <- premiums_ind_ph[idx_updated, ]
P_biv <- premiums_biv_ph[idx_updated, ]

RMSE_ind <- sqrt(mean((P_ind - P_true)^2))
RMSE_biv <- sqrt(mean((P_biv - P_true)^2))

RRMSE_ind <- sqrt(mean(((P_ind - P_true) / P_true)^2))
RRMSE_biv <- sqrt(mean(((P_biv - P_true) / P_true)^2))

round(RRMSE_ind * 100, 2)
round(RRMSE_biv * 100, 2)


## Simulation - rho = -0.25
rho <- -0.25
M <- 2500
set.seed(1)
dep_dat <- sim_dep_gamma(M, rho, 2, 1 / 2, 10 / 3, 10 / 3)

theta1 <- dep_dat$theta1
theta2 <- dep_dat$theta2
n <- dep_dat$n
w <- dep_dat$w
y <- dep_dat$y

round(mean(n == 0), 4)

biv_gam_moments(rho, 2, 1 / 2, 10 / 3, 10 / 3)

M_prem <- 100000
set.seed(1)
prem_theta <- sim_latent_dep_gamma(M_prem, rho, 2, 1 / 2, 10 / 3, 10 / 3)

## Independent fitting
# Frequency
p <- 2
set.seed(123)
mph_ini <- mp(phasetype(structure = "gcoxian", dimension = p))
mph_fit <- fit(mph_ini, n, stepsEM = 10, stepsPH = 50, every = 10, truncationpoint = 20)
#-6155.629

# Severity
wp <- w[n > 0]
np <- n[n > 0]

p <- 3
set.seed(1)
ph_ini <- phasetype(structure = "gcoxian", dimension = p)

ph_fit <- fit_ph(ph_ini, np, wp, stepsEM = 100, stepsPH = 50, every = 1)
#-13487.84

ph_fit@pars$alpha %*% (-ph_fit@pars$S) %*% rep(1, p)


## Joint frequency-severity modeling
p1 <- 6
p2 <- 6
set.seed(19)
bivph_ini <- bivphasetype(dimensions = c(p1, p2))

bmph_fit <- fit_bph(
  bivph_ini,
  n,
  w,
  stepsEM = 100,
  stepsPH = 100,
  nq1 = 200,
  nq2 = 200,
  upper1 = 25,
  upper2 = 8,
  keep_mass = 0.9995,
  keep_ratio = 0.9995,
  every = 1
)

fit_biv <- bmph_fit$bph_fit

bmph_fit$loglik_trace[length(bmph_fit$loglik_trace)]

ph1 <- marginal(fit_biv, 1)
fbph_fit <- mp(ph1)

ph2 <- marginal(fit_biv, 2)
sbph_fit <- frailty(ph2, bhaz = "exponential", bhaz_pars = 1)


# Moments
ph2@pars$alpha %*% (-ph2@pars$S) %*% rep(1, p2)
biv_ph_moments(fit_biv)
round(corr(fit_biv), 4)

# Premiums
premiums_ind_ph <- prem_ph_mat(nclaims, years, wo, mph_fit, ph_fit)
premiums_biv_ph <- prem_bivph_mat(nclaims, years, wo, fit_biv)
premiums_gamma <- prem_gamma_dep_mat(nclaims, years, wo, prem_theta$theta1, prem_theta$theta2)

idx_updated <- 2:nrow(premiums_gamma)

P_true <- premiums_gamma[idx_updated, ]
P_ind <- premiums_ind_ph[idx_updated, ]
P_biv <- premiums_biv_ph[idx_updated, ]

RMSE_ind <- sqrt(mean((P_ind - P_true)^2))
RMSE_biv <- sqrt(mean((P_biv - P_true)^2))

RRMSE_ind <- sqrt(mean(((P_ind - P_true) / P_true)^2))
RRMSE_biv <- sqrt(mean(((P_biv - P_true) / P_true)^2))

round(RRMSE_ind * 100, 2)
round(RRMSE_biv * 100, 2)


## Simulation - rho = 0.5
rho <- -0.5
M <- 2500
set.seed(1)
dep_dat <- sim_dep_gamma(M, rho, 2, 1 / 2, 10 / 3, 10 / 3)

theta1 <- dep_dat$theta1
theta2 <- dep_dat$theta2
n <- dep_dat$n
w <- dep_dat$w
y <- dep_dat$y

round(mean(n == 0), 4)

biv_gam_moments(rho, 2, 1 / 2, 10 / 3, 10 / 3)

M_prem <- 100000
set.seed(1)
prem_theta <- sim_latent_dep_gamma(M_prem, rho, 2, 1 / 2, 10 / 3, 10 / 3)

## Independent fitting
# Frequency
p <- 2
set.seed(123)
mph_ini <- mp(phasetype(structure = "gcoxian", dimension = p))
mph_fit <- fit(mph_ini, n, stepsEM = 10, stepsPH = 50, every = 10, truncationpoint = 20)
#-6155.352

# Severity
wp <- w[n > 0]
np <- n[n > 0]

p <- 3
set.seed(1)
ph_ini <- phasetype(structure = "gcoxian", dimension = p)

ph_fit <- fit_ph(ph_ini, np, wp, stepsEM = 100, stepsPH = 50, every = 1)
#-14552.46
fph_fit <- frailty(ph_fit, bhaz = "exponential", bhaz_pars = 1)


ph_fit@pars$alpha %*% (-ph_fit@pars$S) %*% rep(1, p)


## Joint frequency-severity modeling
p1 <- 6
p2 <- 6
set.seed(1)
bivph_ini <- bivphasetype(dimensions = c(p1, p2))

bmph_fit <- fit_bph(
  bivph_ini,
  n,
  w,
  stepsEM = 50,
  stepsPH = 100,
  nq1 = 200,
  nq2 = 200,
  upper1 = 25,
  upper2 = 8,
  keep_mass = 0.9995,
  keep_ratio = 0.9995,
  every = 1
)

fit_biv <- bmph_fit$bph_fit

bmph_fit$loglik_trace[length(bmph_fit$loglik_trace)]

ph1 <- marginal(fit_biv, 1)
fbph_fit <- mp(ph1)

# Frequency plot
pdf("freq_hist_dep_m05.pdf", width = 6, height = 6)
plot(table(n) / M,
  lwd = 2, ylim = c(0, 0.17),
  ylab = "", xlab = "", main = "Histogram vs fitted density \n (Frequency - Rho = -0.5)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5, col = "darkgray"
)
points(0:max(n), gamma_mix_densf(0:max(n), 2, 0.5), col = "#3498DB", lw = 2, pch = 3, cex = 1.2)
points(0:max(n), dens(mph_fit, 0:max(n)), col = "#339966", lw = 2, cex = 1.2, pch = 2)
points(0:max(n), dens(fbph_fit, 0:max(n)), col = "#E74C3C", lw = 2, cex = 1.2)
legend("topright", inset = 0.05, bty = "n", c("Original density", "Independent PH", "Bivariate PH"), horiz = FALSE, pch = c(3, 2, 1), col = c("#3498DB", "#339966", "#E74C3C"), cex = 1.2)
dev.off()

ph2 <- marginal(fit_biv, 2)
sbph_fit <- frailty(ph2, bhaz = "exponential", bhaz_pars = 1)

# Severity plot
sq <- seq(0, 8, by = 0.01)
pdf("sev_cdf_dep_m05.pdf", width = 6, height = 6)
Fhat <- ecdf(y)
plot(Fhat,
  verticals = TRUE, do.points = FALSE,
  xlab = "", ylab = "",
  main = "Empirical CDF vs fitted CDF \n (Severity - Rho = -0.5)",
  cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5,
  xlim = range(sq), ylim = c(0, 1), lwd = 2
)
lines(sq, gamma_mix_cdf(sq, 10 / 3, 10 / 3), col = "#3498DB", lwd = 2)
lines(sq, cdf(fph_fit, sq), col = "#339966", lwd = 2, lty = 3)
lines(sq, cdf(sbph_fit, sq), col = "#E74C3C", lwd = 2, lty = 5)
legend("bottomright",
  inset = 0.05, bty = "n",
  legend = c("Empirical CDF", "Original CDF", "Independent PH", "Bivariate PH"),
  lty = c(1, 1, 3, 5), lwd = 2,
  col = c("black", "#3498DB", "#339966", "#E74C3C"),
  cex = 1.2
)
dev.off()


# Moments
ph2@pars$alpha %*% (-ph2@pars$S) %*% rep(1, p2)
biv_ph_moments(fit_biv)
round(corr(fit_biv), 4)

# Premiums
premiums_ind_ph <- prem_ph_mat(nclaims, years, wo, mph_fit, ph_fit)
premiums_biv_ph <- prem_bivph_mat(nclaims, years, wo, fit_biv)
premiums_gamma <- prem_gamma_dep_mat(nclaims, years, wo, prem_theta$theta1, prem_theta$theta2)

idx_updated <- 2:nrow(premiums_gamma)

P_true <- premiums_gamma[idx_updated, ]
P_ind <- premiums_ind_ph[idx_updated, ]
P_biv <- premiums_biv_ph[idx_updated, ]

RMSE_ind <- sqrt(mean((P_ind - P_true)^2))
RMSE_biv <- sqrt(mean((P_biv - P_true)^2))

RRMSE_ind <- sqrt(mean(((P_ind - P_true) / P_true)^2))
RRMSE_biv <- sqrt(mean(((P_biv - P_true) / P_true)^2))

round(RRMSE_ind * 100, 2)
round(RRMSE_biv * 100, 2)
