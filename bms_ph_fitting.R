fit_bph_uni <- function(x,
                        n,
                        w,
                        stepsEM = 5,
                        stepsPH = 20,
                        initialpoint1 = 0.0001,
                        truncationpoint1 = 8,
                        delta1 = 0.015,
                        initialpoint2 = 0.0001,
                        truncationpoint2 = 8,
                        delta2 = 0.015,
                        every = 1) {
  start_time <- Sys.time()

  # Biv PH parameters
  x_par <- x@pars
  alpha_fit <- clone_vector(x_par$alpha)
  S11_fit <- clone_matrix(x_par$S11) #* 0.2
  S12_fit <- clone_matrix(x_par$S12) #* 0.2
  S22_fit <- clone_matrix(x_par$S22)

  # Build grid
  grid1 <- seq(initialpoint1, truncationpoint1, by = delta1)
  grid2 <- seq(initialpoint2, truncationpoint2, by = delta2)
  value <- expand.grid(z1 = grid1, z2 = grid2)
  Z <- as.matrix(value)
  Ngrid <- nrow(Z)

  posterior_grid_weights <- function(alpha, S11, S12, S22, n, w, fNY) {
    z1 <- Z[, 1]
    z2 <- Z[, 2]

    # fTheta at grid
    fTheta <- bivph_density(Z, alpha, S11, S12, S22)
    M <- length(n)
    wts <- numeric(Ngrid)

    lgfact <- 2 * lgamma(n + 1L)

    for (g in seq_len(Ngrid)) {
      term <- n * (log(z1[g]) + log(z2[g])) - z1[g] - w * z2[g] - lgfact
      wts[g] <- mean(exp(term) / fNY)
    }
    wts <- wts * fTheta
    wts
  }

  for (k in 1:stepsEM) {
    aux <- bm_cor_dens(n, w, alpha_fit, S11_fit, S12_fit, S22_fit)

    # Discretization of density
    prob <- posterior_grid_weights(alpha_fit, S11_fit, S12_fit, S22_fit, n, w, aux)


    # PH fitting
    for (l in 1:stepsPH) {
      EMstep_bivph_omp(alpha_fit, S11_fit, S12_fit, S22_fit, Z, prob)
    }

    if (k %% every == 0) {
      x_fit <- bivphasetype(alpha = alpha_fit, S11 = S11_fit, S12 = S12_fit, S22 = S22_fit)
      s1 <- sum(lgamma(n + 1))
      s2 <- sum(log(bm_cor_dens(n, w, alpha_fit, S11_fit, S12_fit, S22_fit)))
      cat("\r", "iteration:", k,
        ", logLik:", s1 + s2,
        ", prob disc:", sum(prob * delta1 * delta2),
        ", corr:", corr(x_fit),
        sep = " "
      )
    }
  }


  out_list <- list(bph_fit = x_fit)

  end_time <- Sys.time()
  cat("\n Running time: ", end_time - start_time)

  return(out_list)
}


fit_ph <- function(x,
                   n,
                   w,
                   stepsEM = 5,
                   stepsPH = 20,
                   initialpoint = 0.0001,
                   truncationpoint = 10,
                   maxprobability = 0.01,
                   maxdelta = 0.05,
                   every = 1) {
  start_time <- Sys.time()

  # PH parameters
  x_par <- x@pars
  alpha_fit <- clone_vector(x_par$alpha)
  S_fit <- clone_matrix(x_par$S)

  conditional_density <- function(t, alpha, S, n, w) {
    fTheta <- ph_density(t, alpha, S)
    mean(exp(n * log(t) - w * t - lgamma(n + 1)) * fTheta / ph_laplace_der_vec(w, n + 1, alpha, S))
  }

  for (k in 1:stepsEM) {
    # Discretization of density
    deltat <- 0
    t <- initialpoint

    prob <- numeric(0)
    value <- numeric(0)

    j <- 1

    while (t < truncationpoint) {
      if (conditional_density(t, alpha_fit, S_fit, n, w) < maxprobability / maxdelta) {
        deltat <- maxdelta
      } else {
        deltat <- maxprobability / conditional_density(t, alpha_fit, S_fit, n, w)
      }
      proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, n, w) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, n, w) + conditional_density(t + deltat, alpha_fit, S_fit, n, w))
      while (proba_aux > maxprobability) {
        deltat <- deltat * 0.9
        proba_aux <- deltat / 6 * (conditional_density(t, alpha_fit, S_fit, n, w) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, n, w) + conditional_density(t + deltat, alpha_fit, S_fit, n, w))
      }
      if (proba_aux > 0) {
        value[j] <- (t * conditional_density(t, alpha_fit, S_fit, n, w) + 4 * (t + deltat / 2) * conditional_density(t + deltat / 2, alpha_fit, S_fit, n, w) + (t + deltat) * conditional_density(t + deltat, alpha_fit, S_fit, n, w)) / (conditional_density(t, alpha_fit, S_fit, n, w) + 4 * conditional_density(t + deltat / 2, alpha_fit, S_fit, n, w) + conditional_density(t + deltat, alpha_fit, S_fit, n, w))
        prob[j] <- proba_aux
        j <- j + 1
      }
      t <- t + deltat
    }

    # PH fitting
    for (l in 1:stepsPH) {
      EMstep(alpha_fit, S_fit, value, prob)
      # EMstep_omp(alpha_fit, S_fit, value, prob)
    }

    if (k %% every == 0) {
      s1 <- sum(lgamma(n + 1))
      s2 <- sum(log(ph_laplace_der_vec(w, n + 1, alpha_fit, S_fit)))
      cat("\r", "iteration:", k,
        ", logLik:", s1 + s2,
        ", prob disc:", sum(prob),
        sep = " "
      )
    }
  }

  x_fit <- phasetype(alpha = alpha_fit, S = S_fit)

  end_time <- Sys.time()
  cat("\n Running time: ", end_time - start_time)

  return(x_fit)
}


gauss_legendre <- function(m, a, b) {
  i <- seq_len(m - 1)
  beta <- i / sqrt(4 * i^2 - 1)

  J <- matrix(0, m, m)
  J[cbind(i, i + 1)] <- beta
  J[cbind(i + 1, i)] <- beta

  eig <- eigen(J, symmetric = TRUE)

  nodes <- eig$values
  weights <- 2 * eig$vectors[1, ]^2

  ord <- order(nodes)

  nodes <- nodes[ord]
  weights <- weights[ord]

  list(
    nodes = (b - a) / 2 * nodes + (a + b) / 2,
    weights = (b - a) / 2 * weights
  )
}

make_log_product_quad <- function(nq1 = 80,
                                  nq2 = 80,
                                  lower1 = 1e-5,
                                  upper1 = 20,
                                  lower2 = 1e-5,
                                  upper2 = 6) {
  q1 <- gauss_legendre(nq1, log(lower1), log(upper1))
  q2 <- gauss_legendre(nq2, log(lower2), log(upper2))

  grid1 <- exp(q1$nodes)
  grid2 <- exp(q2$nodes)

  # Jacobian correction for theta = exp(u)
  w1 <- q1$weights * grid1
  w2 <- q2$weights * grid2

  value <- expand.grid(z1 = grid1, z2 = grid2)
  Z <- as.matrix(value)

  quad_weights <- rep(w1, times = length(grid2)) *
    rep(w2, each = length(grid1))

  list(
    Z = Z,
    weight = quad_weights,
    grid1 = grid1,
    grid2 = grid2,
    w1 = w1,
    w2 = w2
  )
}

posterior_quad_weights_grouped <- function(alpha,
                                           S11,
                                           S12,
                                           S22,
                                           n,
                                           w,
                                           A,
                                           Q) {
  Z <- Q$Z
  z1 <- Z[, 1]
  z2 <- Z[, 2]

  M <- length(n)

  fTheta <- bivph_density(Z, alpha, S11, S12, S22)
  fTheta <- pmax(fTheta, 0)

  if (any(A <= 0) || any(!is.finite(A))) {
    stop("bm_cor_dens returned non-positive or non-finite values.")
  }

  base <- numeric(nrow(Z))

  for (r in sort(unique(n))) {
    idx <- which(n == r)

    logc <- -2 * lgamma(r + 1) - log(A[idx]) - log(M)

    # B_r(z2) = sum_{k:n_k=r} exp(-w_k z2) / {M (r!)^2 A_k}
    Br <- vapply(Q$grid2, function(zz) {
      lt <- logc - w[idx] * zz
      m <- max(lt)
      exp(m) * sum(exp(lt - m))
    }, numeric(1))

    Br_rep <- rep(Br, each = length(Q$grid1))

    base <- base +
      exp(-z1 + r * (log(z1) + log(z2))) * Br_rep
  }

  raw_mass <- fTheta * base * Q$weight
  raw_mass <- pmax(raw_mass, 0)

  disc_mass <- sum(raw_mass)

  if (!is.finite(disc_mass) || disc_mass <= 0) {
    stop("Invalid quadrature posterior mass.")
  }

  prob <- raw_mass / disc_mass

  list(
    prob = prob,
    disc_mass = disc_mass
  )
}

prune_quad_nodes <- function(Z,
                             prob,
                             keep_mass = 0.999,
                             keep_ratio = 0.999,
                             eps_theta2 = 1e-6) {
  prob <- pmax(prob, 0)
  prob <- prob / sum(prob)

  take_top <- function(v, target) {
    v <- pmax(v, 0)

    if (sum(v) <= 0) {
      return(integer(0))
    }

    v <- v / sum(v)
    ord <- order(v, decreasing = TRUE)
    cs <- cumsum(v[ord])
    ord[seq_len(which(cs >= target)[1])]
  }

  keep1 <- take_top(prob, keep_mass)

  ratio_contrib <- prob * Z[, 1] / pmax(Z[, 2], eps_theta2)
  keep2 <- take_top(ratio_contrib, keep_ratio)

  keep <- sort(unique(c(keep1, keep2)))

  list(
    Z = Z[keep, , drop = FALSE],
    prob = prob[keep] / sum(prob[keep]),
    keep = keep
  )
}

fit_bph <- function(x,
                    n,
                    w,
                    stepsEM = 100,
                    stepsPH = 50,
                    nq1 = 80,
                    nq2 = 80,
                    lower1 = 1e-5,
                    upper1 = 20,
                    lower2 = 1e-5,
                    upper2 = 6,
                    keep_mass = 0.999,
                    keep_ratio = 0.999,
                    every = 1) {
  start_time <- Sys.time()

  stopifnot(length(n) == length(w))
  stopifnot(all(n >= 0))
  stopifnot(all(w >= 0))
  stopifnot(all(w[n == 0] == 0))

  x_par <- x@pars

  alpha_fit <- clone_vector(x_par$alpha)
  S11_fit <- clone_matrix(x_par$S11)
  S12_fit <- clone_matrix(x_par$S12)
  S22_fit <- clone_matrix(x_par$S22)

  Q <- make_log_product_quad(
    nq1 = nq1,
    nq2 = nq2,
    lower1 = lower1,
    upper1 = upper1,
    lower2 = lower2,
    upper2 = upper2
  )

  loglik_trace <- numeric(stepsEM)
  disc_mass_trace <- numeric(stepsEM)
  n_nodes_trace <- integer(stepsEM)
  corr_trace <- numeric(stepsEM)

  for (k in seq_len(stepsEM)) {
    A <- bm_cor_dens(
      n,
      w,
      alpha_fit,
      S11_fit,
      S12_fit,
      S22_fit
    )

    post <- posterior_quad_weights_grouped(
      alpha_fit,
      S11_fit,
      S12_fit,
      S22_fit,
      n,
      w,
      A,
      Q
    )

    disc_mass_trace[k] <- post$disc_mass

    if (post$disc_mass < 0.995) {
      warning(
        "Quadrature posterior mass below 0.995 at iteration ",
        k,
        ". Consider increasing upper bounds or nq1/nq2."
      )
    }

    pruned <- prune_quad_nodes(
      Z = Q$Z,
      prob = post$prob,
      keep_mass = keep_mass,
      keep_ratio = keep_ratio
    )

    n_nodes_trace[k] <- nrow(pruned$Z)

    for (l in seq_len(stepsPH)) {
      EMstep_bivph_omp(
        alpha_fit,
        S11_fit,
        S12_fit,
        S22_fit,
        pruned$Z,
        pruned$prob
      )
    }

    A_new <- bm_cor_dens(
      n,
      w,
      alpha_fit,
      S11_fit,
      S12_fit,
      S22_fit
    )

    loglik_trace[k] <- sum(lgamma(n + 1) + log(A_new))

    x_fit <- bivphasetype(
      alpha = alpha_fit,
      S11 = S11_fit,
      S12 = S12_fit,
      S22 = S22_fit
    )

    corr_trace[k] <- corr(x_fit)

    if (k %% every == 0) {
      cat(
        "\r",
        "iteration:", k,
        ", logLik:", loglik_trace[k],
        ", quad mass:", disc_mass_trace[k],
        ", nodes:", n_nodes_trace[k],
        ", corr:", corr_trace[k],
        sep = " "
      )
    }
  }

  x_fit <- bivphasetype(
    alpha = alpha_fit,
    S11 = S11_fit,
    S12 = S12_fit,
    S22 = S22_fit
  )

  end_time <- Sys.time()

  cat("\n Running time:", end_time - start_time, "\n")

  list(
    bph_fit = x_fit,
    loglik_trace = loglik_trace,
    disc_mass_trace = disc_mass_trace,
    n_nodes_trace = n_nodes_trace,
    corr_trace = corr_trace
  )
}
