### Independent PH components ###

inverse_ph_moment <- function(ph_mod) {
  s <- -ph_mod@pars$S %*% rep(1, p)
  aux <- -ph_mod@pars$alpha %*% logm(-ph_mod@pars$S) %*% s
  aux[1, 1]
}


# Premiums with independent model
prem_ph <- function(n, E, w, alpha1, S1, alpha2, S2) {
  p1 <- length(alpha1)
  p2 <- length(alpha2)
  I1 <- diag(p1)
  I2 <- diag(p2)
  s1 <- -S1 %*% rep(1, p1)
  s2 <- -S2 %*% rep(1, p2)
  inv_S1 <- solve(E * I1 - S1)
  inv_S2 <- solve(w * I2 - S2)
  prem <- 0
  if (n == 0 && E == 0) {
    aux_freq1 <- alpha1 %*% inv_S1 %*% rep(1, p1)
    aux_sev1 <- -alpha2 %*% logm(-S2) %*% s2
    prem <- (aux_freq1[1, 1]) * aux_sev1[1, 1]
  } else if (n == 0 && E > 0) {
    inv_pow1 <- matrix_power(n + 2, inv_S1)
    aux_freq1 <- alpha1 %*% inv_pow1 %*% s1
    aux_freq2 <- alpha1 %*% inv_pow1 %*% (E * I1 - S1) %*% s1
    aux_sev1 <- -alpha2 %*% logm(-S2) %*% s2
    prem <- (aux_freq1[1, 1] / aux_freq2[1, 1]) * aux_sev1[1, 1]
  } else {
    inv_pow1 <- matrix_power(n + 2, inv_S1)
    inv_pow2 <- matrix_power(n + 1, inv_S2)
    aux_freq1 <- alpha1 %*% inv_pow1 %*% s1
    aux_freq2 <- alpha1 %*% inv_pow1 %*% (E * I1 - S1) %*% s1
    aux_sev1 <- alpha2 %*% inv_pow2 %*% s2
    aux_sev2 <- alpha2 %*% inv_pow2 %*% (w * I2 - S2) %*% s2
    prem <- ((n + 1) / n) * (aux_freq1[1, 1] / aux_freq2[1, 1]) * (aux_sev2[1, 1] / aux_sev1[1, 1])
  }
  prem
}

prem_ph_mat <- function(nclaims, years, w, ph_freq, ph_sev) {
  alpha1 <- ph_freq@pars$alpha
  S1 <- ph_freq@pars$S
  alpha2 <- ph_sev@pars$alpha
  S2 <- ph_sev@pars$S

  premiums <- matrix(0, length(years) + 1, length(nclaims) + 1)
  for (t in c(0, years)) {
    premiums[t + 1, 1] <- prem_ph(0, t, 0, alpha1, S1, alpha2, S2)
  }

  for (t in years) {
    for (ni in nclaims) {
      premiums[t + 1, ni + 1] <- prem_ph(ni, t, wo, alpha1, S1, alpha2, S2)
    }
  }
  premiums
}


### Bivariate PH ###

# Moments
biv_ph_moments <- function(biv_ph) {
  ph1 <- marginal(biv_ph, 1)
  ph2 <- marginal(biv_ph, 2)

  m1 <- moment(ph1, 1)

  s <- -ph2@pars$S %*% rep(1, p2)
  m2 <- -ph2@pars$alpha %*% logm(-ph2@pars$S) %*% s

  eta <- biv_ph@pars$alpha
  S11 <- biv_ph@pars$S11
  S12 <- biv_ph@pars$S12
  S22 <- biv_ph@pars$S22
  s2 <- -S22 %*% rep(1, p2)
  m12 <- -eta %*% matrix_power(2, solve(-S11)) %*% S12 %*% logm(-S22) %*% s2
  list(m1 = m1, m2 = m2[1, 1], m12 = m12[1, 1], covNY = m12[1, 1] - m1 * m2[1, 1])
}

# Premiums
prem_bivph <- function(n, E, w, eta, S11, S12, S22) {
  p1 <- length(eta)
  p2 <- nrow(S22)
  I1 <- diag(p1)
  I2 <- diag(p2)
  s1 <- -S11 %*% rep(1, p1)
  s2 <- -S22 %*% rep(1, p2)
  inv_S1 <- solve(E * I1 - S11)
  inv_S2 <- solve(w * I2 - S22)
  prem <- 0
  if (n == 0 && E == 0) {
    aux <- -eta %*% matrix_power(2, inv_S1) %*% S12 %*% logm(-S22) %*% s2
    prem <- aux[1, 1]
  } else if (n == 0 && E > 0) {
    inv_pow1 <- matrix_power(n + 2, inv_S1)
    aux_freq1 <- eta %*% inv_pow1 %*% S12 %*% logm(-S22) %*% s2
    aux_freq2 <- eta %*% inv_pow1 %*% (E * I1 - S11) %*% s1
    prem <- -(aux_freq1[1, 1] / aux_freq2[1, 1])
  } else {
    inv_pow1 <- matrix_power(n + 2, inv_S1)
    inv_pow2 <- matrix_power(n + 1, inv_S2)
    aux_freq <- eta %*% inv_pow1 %*% S12 %*% inv_pow2 %*% (w * I2 - S22) %*% s2
    aux_sev <- eta %*% inv_pow1 %*% (E * I1 - S11) %*% S12 %*% inv_pow2 %*% s2
    prem <- ((n + 1) / n) * (aux_freq[1, 1] / aux_sev[1, 1])
  }
  prem
}

prem_bivph_mat <- function(nclaims, years, w, biv_ph) {
  eta <- biv_ph@pars$alpha
  S11 <- biv_ph@pars$S11
  S12 <- biv_ph@pars$S12
  S22 <- biv_ph@pars$S22

  premiums <- matrix(0, length(years) + 1, length(nclaims) + 1)
  for (t in c(0, years)) {
    premiums[t + 1, 1] <- prem_bivph(0, t, 0, eta, S11, S12, S22)
  }

  for (t in years) {
    for (ni in nclaims) {
      premiums[t + 1, ni + 1] <- prem_bivph(ni, t, wo, eta, S11, S12, S22)
    }
  }
  premiums
}
