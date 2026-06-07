### Independent case ###

# Simulation
sim_ind_gamma <- function(M, rate1, shape1, rate2, shape2) {
  theta1 <- rgamma(M, rate1, shape1)
  theta2 <- rgamma(M, rate2, shape2)
  n <- rep(0, M)
  w <- rep(0, M)
  y <- numeric(0)
  for (i in 1:M) {
    n[i] <- rpois(1, theta1[i])
    if (n[i] > 0) {
      sim_sev <- rexp(n[i], theta2[i])
      y <- c(y, sim_sev)
      w[i] <- sum(sim_sev)
    } else {
      w[i] <- 0
    }
  }
  list(theta1 = theta1, theta2 = theta2, n = n, w = w, y = y)
}

## True densities/cdfs
# Frequency
gamma_mix_densf <- function(n, kappa, lambda) {
  gamma(n + kappa) / (gamma(kappa) * gamma(n + 1)) * (lambda / (lambda + 1))^kappa * (1 / (lambda + 1))^n
}
# Severity
gamma_mix_denss <- function(y, kappa, lambda) {
  (gamma(kappa + 1) / gamma(kappa)) * lambda^kappa * (lambda + y)^(-kappa - 1)
}
gamma_mix_cdf <- function(y, kappa, lambda) {
  1 - lambda^kappa * (lambda + y)^(-kappa)
}
gamma_mix_denssl <- function(n, w, kappa, lambda) {
  (gamma(kappa + n) / gamma(kappa)) * lambda^kappa * (lambda + w)^(-kappa - n)
}

## Premiums
prem_gamma <- function(n, E, w, kappa1, lambda1, kappa2, lambda2) {
  ((kappa1 + n) / (lambda1 + E)) * (lambda2 + w) / (kappa2 + n - 1)
}

prem_gamma_mat <- function(nclaims, years, w, kappa1, lambda1, kappa2, lambda2) {
  premiums <- matrix(0, length(years) + 1, length(nclaims) + 1)
  for (t in c(0, years)) {
    premiums[t + 1, 1] <- prem_gamma(0, t, 0, kappa1, lambda1, kappa2, lambda2)
  }

  for (t in years) {
    for (ni in nclaims) {
      premiums[t + 1, ni + 1] <- prem_gamma(ni, t, w, kappa1, lambda1, kappa2, lambda2)
    }
  }
  premiums
}

# Loglikelihood using ind Gamma
ind_gamma <- function(n, z, kappa1, lambda1, kappa2, lambda2) {
  ll <- rep(0, length(n))
  for (i in 1:length(n)) {
    if (n[i] > 0) {
      aux1 <- log(gamma(n[i] + kappa1) / (gamma(kappa1) * gamma(n[i] + 1)) * (lambda1 / (lambda1 + 1))^kappa1 * (1 / (lambda1 + 1))^n[i])
      aux2 <- log((gamma(kappa2 + n[i]) / gamma(kappa2)) * lambda2^kappa2 * (lambda2 + z[i])^(-kappa2 - n[i]))
      ll[i] <- aux1 + aux2
    } else {
      ll[i] <- log(gamma(n[i] + kappa1) / (gamma(kappa1) * gamma(n[i] + 1)) * (lambda1 / (lambda1 + 1))^kappa1 * (1 / (lambda1 + 1))^n[i])
    }
  }
  sum(ll)
}

### Dependent case ###

# Simulation
sim_dep_gamma <- function(M, rho, shape1, rate1, shape2, rate2) {
  cc <- normalCopula(param = rho, dim = 2)
  u <- rCopula(M, copula = cc)
  theta1 <- qgamma(u[, 1], shape1, rate1)
  theta2 <- qgamma(u[, 2], shape2, rate2)

  n <- rep(0, M)
  w <- rep(0, M)
  y <- numeric(0)
  for (i in 1:M) {
    n[i] <- rpois(1, theta1[i])
    if (n[i] > 0) {
      sim_sev <- rexp(n[i], theta2[i])
      y <- c(y, sim_sev)
      w[i] <- sum(sim_sev)
    } else {
      w[i] <- 0
    }
  }
  list(theta1 = theta1, theta2 = theta2, n = n, w = w, y = y)
}

sim_latent_dep_gamma <- function(M, rho, shape1, rate1, shape2, rate2) {
  cc <- normalCopula(param = rho, dim = 2)
  u <- rCopula(M, copula = cc)
  list(
    theta1 = qgamma(u[, 1], shape = shape1, rate = rate1),
    theta2 = qgamma(u[, 2], shape = shape2, rate = rate2)
  )
}

# Moments
biv_gam_moments <- function(rho, shape1, rate1, shape2, rate2, lb = c(0, 0.0001), ub = c(20, 20)) {
  cc <- normalCopula(param = rho, dim = 2)

  mvm <- mvdc(
    copula = cc,
    margins = c("gamma", "gamma"),
    paramMargins = list(list(shape = shape1, rate = rate1), list(shape = shape2, rate = rate2))
  )

  aux_fn <- function(x) {
    dMvdc(c(x[1], x[2]), mvm) * x[1] / x[2]
  }
  m12 <- adaptIntegrate(aux_fn, lowerLimit = lb, upperLimit = ub)$integral

  m1 <- integrate(function(x) x * dgamma(x, shape1, rate1), 0, Inf)$value
  m2 <- integrate(function(x) dgamma(x, shape2, rate2) / x, 0, Inf)$value
  list(m1 = m1, m2 = m2, m12 = m12, covNY = m12 - m1 * m2)
}

# Premiums using MC
prem_gamma_copula_mc <- function(n, E, w, theta1, theta2, eps = 1e-300) {
  log_weight <- n * (log(theta1) + log(theta2)) - E * theta1 - w * theta2
  m <- max(log_weight)
  weight <- exp(log_weight - m)
  sum((theta1 / theta2) * weight) / pmax(sum(weight), eps)
}

prem_gamma_dep_mat <- function(nclaims, years, w, theta1, theta2, eps = 1e-300) {
  premiums <- matrix(0, length(years) + 1, length(nclaims) + 1)
  for (t in c(0, years)) {
    premiums[t + 1, 1] <- prem_gamma_copula_mc(0, t, 0, theta1, theta2, eps)
  }

  for (t in years) {
    for (ni in nclaims) {
      premiums[t + 1, ni + 1] <- prem_gamma_copula_mc(ni, t, w, theta1, theta2, eps)
    }
  }
  premiums
}
