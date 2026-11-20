# =============================================================================
# Bayesian Inference for Survival Analysis in Two-Threshold Generalized
# δ-shock Models
#
# Companion code (simulation study) for:
#   L. Rajabi and S. Khazaei, "Bayesian Inference for Survival Analysis in
#   Two-Threshold Generalized δ-shock Models",
#   Communications in Statistics - Theory and Methods (submitted).
#
# Reproduces Tables 1-4 and Figures 1-4 of the simulation study.
#
# R >= 4.0, base R only.
# Run from the repository root:  source("01_simulation.R")
# =============================================================================

set.seed(2026)
OUT_DIR <- "outputs_simulation"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------- settings --
p   <- 0.3
d1  <- 1
d2  <- 3
n   <- 150
THETAS <- c(0.1, 0.3, 0.5, 0.7, 0.9)
TS     <- c(5, 10, 20)

R      <- 500
Rc     <- 50
R3     <- 50
R4     <- 50
Bboot  <- 200
Rh     <- 100
NBIG   <- 200
RB     <- 20

q0   <- 1 - p
p2_0 <- q0^d1 - q0^d2
p3_0 <- 1 - q0^d1
mean_T0 <- function(th) (1 + p2_0 * th + p3_0) / (p * (p2_0 * th + p3_0))
TMAX <- max(400, ceiling(25 * mean_T0(0)) + 50)

# ------------------------------------------------------------ model setup --
setup_C <- function(p, d1, d2, tmax) {
  q <- 1 - p
  x <- 1:tmax
  pmf <- p * q^(x - 1)
  det  <- ifelse(x <= d1, pmf, 0)
  prob <- ifelse(x > d1 & x <= d2, pmf, 0)
  non  <- ifelse(x > d2, pmf, 0)
  M <- tmax %/% (d1 + 1) + 1
  g <- C1 <- C2 <- matrix(0, tmax, M + 1)
  g[, 1] <- pmf
  for (s in 1:(tmax - 1)) {
    for (dx in 1:(tmax - s)) {
      sp <- s + dx
      if (non[dx] > 0) g[sp, ] <- g[sp, ] + g[s, ] * non[dx]
      if (prob[dx] > 0) {
        g[sp, 2:(M + 1)] <- g[sp, 2:(M + 1)] + g[s, 1:M] * prob[dx]
        C1[sp, ] <- C1[sp, ] + g[s, ] * prob[dx]
      }
      if (det[dx] > 0) C2[sp, ] <- C2[sp, ] + g[s, ] * det[dx]
    }
  }
  list(C1 = C1, C2 = C2)
}
CC <- setup_C(p, d1, d2, TMAX)

for (th in c(0.1, 0.5, 0.9)) {
  m <- 0:(ncol(CC$C1) - 1)
  f <- sapply(1:TMAX, function(t) sum((1 - th)^m * (th * CC$C1[t, ] + CC$C2[t, ])))
  stopifnot(abs(sum(f) - 1) < 1e-6, abs(sum((1:TMAX) * f) - mean_T0(th)) < 1e-3)
}
cat("pmf check passed\n")

# --------------------------------------------------------- likelihood/post --
Dcache <- vector("list", TMAX)
getD <- function(t) {
  if (is.null(Dcache[[t]])) {
    mm <- max(which(CC$C1[t, ] > 0 | CC$C2[t, ] > 0))
    Dcache[[t]] <<- rbind(CC$C2[t, 1:mm], CC$C1[t, 1:mm])
  }
  Dcache[[t]]
}
conv_rows <- function(cur, y) {
  c0 <- ncol(cur)
  out <- matrix(0, nrow(cur), c0 + length(y) - 1)
  for (j in seq_along(y)) {
    if (y[j] > 0) {
      idx <- j:(j + c0 - 1)
      out[, idx] <- out[, idx] + cur * y[j]
    }
  }
  out
}
coef_array <- function(t) {
  cur <- matrix(1, 1, 1)
  for (ti in t) {
    D <- getD(ti)
    nr <- nrow(cur)
    new <- matrix(0, nr + 1, ncol(cur) + ncol(D) - 1)
    new[1:nr, ] <- new[1:nr, ] + conv_rows(cur, D[1, ])
    new[2:(nr + 1), ] <- new[2:(nr + 1), ] + conv_rows(cur, D[2, ])
    cur <- new / max(new)
  }
  cur
}
make_post <- function(cur) {
  A <- row(cur) - 1
  Mm <- col(cur) - 1
  keep <- cur > 0
  list(lc = log(cur[keep]), a = A[keep], m = Mm[keep])
}
loglik <- function(P, th) {
  sapply(th, function(x) {
    z <- P$lc + P$a * log(x) + P$m * log1p(-x)
    mx <- max(z)
    mx + log(sum(exp(z - mx)))
  })
}
post_w <- function(P, al, be) {
  lw <- P$lc + lbeta(al + P$a, be + P$m)
  w <- exp(lw - max(lw))
  w / sum(w)
}
post_mean <- function(P, al = 1, be = 1) {
  w <- post_w(P, al, be)
  sum(w * (al + P$a) / (al + be + P$a + P$m))
}
post_ci <- function(P, al = 1, be = 1, level = 0.95) {
  w <- post_w(P, al, be)
  k <- w > 1e-12
  w <- w[k] / sum(w[k])
  sh1 <- al + P$a[k]
  sh2 <- be + P$m[k]
  cdf <- function(x) sum(w * pbeta(x, sh1, sh2))
  qf <- function(pr) uniroot(function(x) cdf(x) - pr, c(1e-10, 1 - 1e-10), tol = 1e-8)$root
  c(qf((1 - level) / 2), qf(1 - (1 - level) / 2))
}
mle <- function(P) {
  g <- seq(0.001, 0.999, length.out = 300)
  l <- loglik(P, g)
  i <- which.max(l)
  optimize(function(x) loglik(P, x), c(g[max(i - 1, 1)], g[min(i + 1, 300)]),
           maximum = TRUE, tol = 1e-8)$maximum
}
lr_ci <- function(P, th_hat) {
  g <- seq(0.001, 0.999, length.out = 400)
  l <- loglik(P, g)
  lmax <- loglik(P, th_hat)
  ok <- g[2 * (lmax - l) <= qchisq(0.95, 1)]
  if (length(ok) == 0) return(c(0.001, 0.999))
  c(min(ok), max(ok))
}
est <- function(t, al = 1, be = 1, ci = TRUE, lr = TRUE) {
  P <- make_post(coef_array(t))
  th <- mle(P)
  out <- list(mle = th, bayes = post_mean(P, al, be), P = P)
  if (ci) {
    out$cr <- post_ci(P, al, be)
    if (lr) out$lr <- lr_ci(P, th)
  }
  out
}

# ----------------------------------------------------------- survival func --
make_S <- function(t0) {
  c1 <- colSums(CC$C1[1:t0, , drop = FALSE])
  c2 <- colSums(CC$C2[1:t0, , drop = FALSE])
  m <- 0:(length(c1) - 1)
  function(th) sapply(th, function(x) 1 - sum((1 - x)^m * (x * c1 + c2)))
}
post_S_mean <- function(P, t0, al = 1, be = 1) {
  c1 <- colSums(CC$C1[1:t0, , drop = FALSE])
  c2 <- colSums(CC$C2[1:t0, , drop = FALSE])
  kk <- which(c1 > 0 | c2 > 0) - 1
  w <- post_w(P, al, be)
  keep <- w > 1e-12
  w <- w[keep] / sum(w[keep])
  A <- al + P$a[keep]
  Bb <- be + P$m[keep]
  e <- numeric(length(w))
  for (i in seq_along(w)) {
    e[i] <- sum(c1[kk + 1] * exp(lbeta(A[i] + 1, Bb[i] + kk) - lbeta(A[i], Bb[i])) +
                c2[kk + 1] * exp(lbeta(A[i], Bb[i] + kk) - lbeta(A[i], Bb[i])))
  }
  1 - sum(w * e)
}

# ------------------------------------------------------------- simulation --
sim_lifetime <- function(theta) {
  t <- rgeom(1, p) + 1
  repeat {
    x <- rgeom(1, p) + 1
    t <- t + x
    if (x <= d1) break
    if (x <= d2 && runif(1) < theta) break
  }
  t
}
sim_sample <- function(theta, n) {
  repeat {
    t <- replicate(n, sim_lifetime(theta))
    if (all(t <= TMAX)) return(t)
  }
}
sim_het <- function(n, theta) {
  repeat {
    t <- replicate(n, {
      pp <- sample(c(0.15, 0.6), 1)
      tt <- rgeom(1, pp) + 1
      repeat {
        x <- rgeom(1, pp) + 1
        tt <- tt + x
        if (x <= d1) break
        if (x <= d2 && runif(1) < theta) break
      }
      tt
    })
    if (all(t <= TMAX)) return(t)
  }
}

# --------------------------------------------------- predictive check (HPC) --
pmf_theta <- function(x) {
  m <- 0:(ncol(CC$C1) - 1)
  as.vector((x * CC$C1 + CC$C2) %*% (1 - x)^m)
}
sample_theta <- function(P, k, al = 1, be = 1) {
  w <- post_w(P, al, be)
  i <- sample.int(length(w), k, TRUE, w)
  rbeta(k, al + P$a[i], be + P$m[i])
}
hpc_p2 <- function(t, q = 0.5, nrep = Rh) {
  N <- length(t)
  nobs <- ceiling(q * N)
  idx <- sample.int(N)
  obs <- t[idx[1:nobs]]
  ho <- t[idx[(nobs + 1):N]]
  th <- sample_theta(make_post(coef_array(obs)), nrep)
  cnt <- 0
  for (x in th) {
    pm <- pmf_theta(x)
    yrep <- sample.int(TMAX, length(ho), TRUE, pm)
    d_rep <- sd(yrep)
    d_new <- sd(ho)
    cnt <- cnt + (d_rep > d_new) + 0.5 * (d_rep == d_new)
  }
  pp <- cnt / nrep
  2 * min(pp, 1 - pp)
}

# ------------------------------------------------------ bootstrap intervals --
grid_b <- seq(0.001, 0.999, length.out = 999)
obs_loglik <- function(t) {
  W <- outer(1 - grid_b, 0:(ncol(CC$C1) - 1), "^")
  A <- W %*% t(CC$C1[t, , drop = FALSE])
  B0 <- W %*% t(CC$C2[t, , drop = FALSE])
  log(pmax(grid_b * A + B0, 1e-300))
}
boot_ci <- function(LL, type) {
  nn <- ncol(LL)
  if (type == "np") {
    Wt <- rmultinom(Bboot, nn, rep(1 / nn, nn))
  } else {
    x <- matrix(rexp(nn * Bboot), nn, Bboot)
    Wt <- sweep(x, 2, colSums(x), "/") * nn
  }
  est_b <- grid_b[apply(LL %*% Wt, 2, which.max)]
  as.numeric(quantile(est_b, c(0.025, 0.975)))
}

# ---------------------------------------------------------- CRLB and Fisher --
fisher_crlb <- function(thetas, n) {
  h <- 1e-4
  out <- data.frame()
  for (th in thetas) {
    P0 <- pmf_theta(th)
    dP <- (pmf_theta(th + h) - pmf_theta(th - h)) / (2 * h)
    k <- P0 > 1e-14
    I <- sum(dP[k]^2 / P0[k])
    out <- rbind(out, data.frame(theta = th, Fisher_info = I,
                                 CRLB_MSE = 1 / (n * I), CRLB_SD = 1 / sqrt(n * I)))
  }
  out
}

# =============================================================================
# Table 1
# =============================================================================
cat("\nTable 1\n")
tab1 <- data.frame()
for (th in THETAS) {
  m <- t(replicate(R, { e <- est(sim_sample(th, n), ci = FALSE); c(e$mle, e$bayes) }))
  mb <- t(replicate(RB, { e <- est(sim_sample(th, NBIG), ci = FALSE); c(e$mle, e$bayes) }))
  cat("  theta =", th, "\n")
  tab1 <- rbind(tab1, data.frame(
    theta = th,
    mean_Bayes = mean(m[, 2]), mean_MLE = mean(m[, 1]),
    Bayes_bias = mean(m[, 2]) - th, Bayes_MSE = mean((m[, 2] - th)^2),
    MLE_bias = mean(m[, 1]) - th, MLE_MSE = mean((m[, 1] - th)^2),
    MCSE_MSE_Bayes = sd((m[, 2] - th)^2) / sqrt(R),
    MCSE_MSE_MLE = sd((m[, 1] - th)^2) / sqrt(R),
    Bayes_MSE_nBIG = mean((mb[, 2] - th)^2),
    MLE_MSE_nBIG = mean((mb[, 1] - th)^2)))
  write.csv(tab1, file.path(OUT_DIR, "table1_partial.csv"), row.names = FALSE)
  flush.console()
}
names(tab1)[names(tab1) == "Bayes_MSE_nBIG"] <- paste0("Bayes_MSE_n", NBIG)
names(tab1)[names(tab1) == "MLE_MSE_nBIG"] <- paste0("MLE_MSE_n", NBIG)
tab1$CRLB_MSE <- fisher_crlb(THETAS, n)$CRLB_MSE

# =============================================================================
# Table 2
# =============================================================================
cat("\nTable 2\n")
tab2 <- data.frame()
for (th in THETAS) {
  cover <- numeric(4)
  width <- matrix(0, Rc, 4)
  rej_ok <- 0
  rej_mis <- 0
  for (r in 1:Rc) {
    tt <- sim_sample(th, n)
    e <- est(tt)
    LL <- obs_loglik(tt)
    ints <- list(e$cr, e$lr, boot_ci(LL, "np"), boot_ci(LL, "wlb"))
    for (k in 1:4) {
      cover[k] <- cover[k] + (ints[[k]][1] <= th & th <= ints[[k]][2])
      width[r, k] <- diff(ints[[k]])
    }
    rej_ok <- rej_ok + (hpc_p2(tt) < 0.05)
    rej_mis <- rej_mis + (hpc_p2(sim_het(n, th)) < 0.05)
  }
  cat("  theta =", th, "\n")
  tab2 <- rbind(tab2, data.frame(
    theta = th,
    cov_credible = cover[1] / Rc, width_credible = mean(width[, 1]),
    cov_LR = cover[2] / Rc, width_LR = mean(width[, 2]),
    cov_NPboot = cover[3] / Rc, width_NPboot = mean(width[, 3]),
    cov_WLboot = cover[4] / Rc, width_WLboot = mean(width[, 4]),
    hpc_reject_correct = rej_ok / Rc, hpc_reject_misspec = rej_mis / Rc))
}

# =============================================================================
# Table 3
# =============================================================================
cat("\nTable 3\n")
priors <- list("Beta(1,1)" = c(1, 1), "Beta(.5,.5)" = c(.5, .5),
               "Beta(2,2)" = c(2, 2), "Beta(5,5)" = c(5, 5),
               "mean .2, ESS 10" = c(2, 8), "mean .8, ESS 10" = c(8, 2))
tab3 <- matrix(NA, length(priors) + 1, length(THETAS),
               dimnames = list(c(names(priors), "MLE"), THETAS))
for (j in seq_along(THETAS)) {
  th <- THETAS[j]
  err <- matrix(NA, R3, length(priors) + 1)
  for (r in 1:R3) {
    P <- make_post(coef_array(sim_sample(th, n)))
    for (k in seq_along(priors)) {
      err[r, k] <- (post_mean(P, priors[[k]][1], priors[[k]][2]) - th)^2
    }
    err[r, length(priors) + 1] <- (mle(P) - th)^2
  }
  tab3[, j] <- colMeans(err)
  cat("  theta =", th, "\n")
}

# =============================================================================
# Table 4
# =============================================================================
cat("\nTable 4\n")
Sfun <- lapply(TS, make_S)
tab4 <- data.frame()
for (th in THETAS) {
  S_true <- sapply(Sfun, function(f) f(th))
  E <- array(NA, c(R4, length(TS), 3))
  CVc <- CVl <- Lc <- Uc <- Ll <- Ul <- matrix(NA, R4, length(TS))
  for (r in 1:R4) {
    e <- est(sim_sample(th, n))
    for (k in seq_along(TS)) {
      E[r, k, 1] <- Sfun[[k]](e$bayes)
      E[r, k, 2] <- Sfun[[k]](e$mle)
      E[r, k, 3] <- post_S_mean(e$P, TS[k])
      ci <- sort(Sfun[[k]](e$cr))
      li <- sort(Sfun[[k]](e$lr))
      CVc[r, k] <- (ci[1] <= S_true[k] & S_true[k] <= ci[2])
      CVl[r, k] <- (li[1] <= S_true[k] & S_true[k] <= li[2])
      Lc[r, k] <- ci[1]; Uc[r, k] <- ci[2]
      Ll[r, k] <- li[1]; Ul[r, k] <- li[2]
    }
  }
  cat("  theta =", th, "\n")
  for (k in seq_along(TS)) {
    err <- E[, k, ] - S_true[k]
    tab4 <- rbind(tab4, data.frame(
      theta = th, t = TS[k], S_true = S_true[k],
      S_Bayes_plugin = mean(E[, k, 1]), S_MLE_plugin = mean(E[, k, 2]), S_Bayes_post = mean(E[, k, 3]),
      MSE_Bayes_plugin = mean(err[, 1]^2), MSE_MLE_plugin = mean(err[, 2]^2), MSE_Bayes_post = mean(err[, 3]^2),
      cov_credible = mean(CVc[, k]), cov_LR = mean(CVl[, k]),
      cred_lower = mean(Lc[, k]), cred_upper = mean(Uc[, k]),
      LR_lower = mean(Ll[, k]), LR_upper = mean(Ul[, k])))
  }
}
tab4_paper <- tab4[, c("theta", "t", "S_true", "S_Bayes_post", "S_MLE_plugin",
                       "MSE_Bayes_post", "MSE_MLE_plugin",
                       "cred_lower", "cred_upper", "LR_lower", "LR_upper",
                       "cov_credible", "cov_LR")]
names(tab4_paper)[names(tab4_paper) == "S_Bayes_post"] <- "S_Bayes"
names(tab4_paper)[names(tab4_paper) == "S_MLE_plugin"] <- "S_MLE"
names(tab4_paper)[names(tab4_paper) == "MSE_Bayes_post"] <- "MSE_Bayes"
names(tab4_paper)[names(tab4_paper) == "MSE_MLE_plugin"] <- "MSE_MLE"

# =============================================================================
# Figures
# =============================================================================
RUN_ID <- format(Sys.time(), "%H%M%S")
show_and_save <- function(name, w, h, draw) {
  if (interactive()) {
    dev.new(width = w / 100, height = h / 100)
    draw()
  }
  png(file.path(OUT_DIR, paste0(name, "_", RUN_ID, ".png")), w * 2, h * 2, res = 200)
  draw()
  dev.off()
}

th_fig <- 0.5
set.seed(2026 + 1)
R_fig <- 1000
m5 <- t(replicate(R_fig, { e <- est(sim_sample(th_fig, n), ci = FALSE); c(e$mle, e$bayes) }))
Ks <- c(500, 1000)

fig1 <- function() {
  par(mfrow = c(length(Ks), 2), mar = c(4, 4, 3, 1), cex.main = 0.9)
  for (K in Ks) {
    yl <- range(m5[1:K, 2], m5[1:K, 1], th_fig)
    plot(m5[1:K, 1], type = "l", col = "purple", ylim = yl,
         xlab = "replication", ylab = expression(theta),
         main = bquote("MLE" ~ theta ~ "estimates" ~ (R == .(K))))
    abline(h = th_fig, col = "red", lwd = 2)
    plot(m5[1:K, 2], type = "l", col = "blue", ylim = yl,
         xlab = "replication", ylab = expression(theta),
         main = bquote("Bayesian" ~ theta ~ "estimates" ~ (R == .(K))))
    abline(h = th_fig, col = "red", lwd = 2)
  }
}
fig2 <- function() {
  par(mfrow = c(length(Ks), 2), mar = c(4, 4, 3, 1), cex.main = 0.9)
  for (K in Ks) {
    for (j in 1:2) {
      rm <- cumsum(m5[1:K, j]) / (1:K)
      yl <- range(c(rm, th_fig, mean(m5[1:K, j])))
      yl <- yl + c(-0.02, 0.02) * diff(yl)
      plot(rm, type = "l", col = c("purple", "blue")[j], lwd = 2, ylim = yl,
           yaxt = "n", xlab = "replication", ylab = "running mean",
           main = paste0("Running mean, ", c("MLE", "Bayes")[j], " (R = ", K, ")"))
      axis(2, las = 1)
      abline(h = th_fig, col = "red", lwd = 2, lty = 2)
      abline(h = mean(m5[1:K, j]), col = "darkgreen", lwd = 2, lty = 2)
      legend("topright", c("running mean", "true theta", "mean of estimates"),
             col = c(c("purple", "blue")[j], "red", "darkgreen"),
             lty = c(1, 2, 2), lwd = 2, bty = "n", cex = 0.8)
    }
  }
}

mean_T_p <- function(pp, d1, d2, theta) {
  qq <- 1 - pp
  pp2 <- qq^d1 - qq^d2
  pp3 <- 1 - qq^d1
  (1 + pp2 * theta + pp3) / (pp * (pp2 * theta + pp3))
}
p_grid <- seq(0.1, 0.9, by = 0.05)
theta_fixed <- 0.5
R_f3 <- 50
et_p <- matrix(NA, length(p_grid), 3, dimnames = list(NULL, c("True", "Bayes", "MLE")))

p_save <- p; TMAX_save <- TMAX; CC_save <- CC; Dcache_save <- Dcache
for (i in seq_along(p_grid)) {
  p <- p_grid[i]
  TMAX <- max(400, ceiling(15 * mean_T_p(p, d1, d2, 0.9)) + 50)
  CC <- setup_C(p, d1, d2, TMAX)
  Dcache <- vector("list", TMAX)
  eB <- numeric(R_f3)
  eM <- numeric(R_f3)
  for (r in 1:R_f3) {
    P <- make_post(coef_array(sim_sample(theta_fixed, n)))
    eB[r] <- post_mean(P)
    eM[r] <- mle(P)
  }
  et_p[i, ] <- c(mean_T_p(p, d1, d2, theta_fixed),
                 mean(mean_T_p(p, d1, d2, eB)),
                 mean(mean_T_p(p, d1, d2, eM)))
  cat("  E(T) figure: p =", p, "\n")
}
p <- p_save; TMAX <- TMAX_save; CC <- CC_save; Dcache <- Dcache_save

fig3 <- function() {
  par(mar = c(4.5, 4.5, 3, 1), cex.main = 0.9)
  plot(p_grid, et_p[, "True"], type = "l", lwd = 2, col = "red",
       xlab = "p", ylab = "E(T)",
       main = bquote(E(T) ~ "vs" ~ p ~ (delta[1] == .(d1) * "," ~ delta[2] == .(d2) * "," ~ theta == .(theta_fixed))),
       ylim = range(et_p, na.rm = TRUE))
  lines(p_grid, et_p[, "Bayes"], lty = 2, col = "blue", lwd = 2)
  lines(p_grid, et_p[, "MLE"], lty = 3, col = "darkgreen", lwd = 2)
  legend("topright", c("True", "Bayes estimate", "MLE estimate"),
         col = c("red", "blue", "darkgreen"), lty = c(1, 2, 3), lwd = 2, bty = "n")
}
fig_surv <- function() {
  t_seq <- seq(1, min(60, TMAX), by = 1)
  S_vec <- function(th) sapply(t_seq, function(t0) make_S(t0)(th))
  S_true <- S_vec(th_fig)
  S_bayes <- S_vec(mean(m5[1:R, 2]))
  S_mle   <- S_vec(mean(m5[1:R, 1]))
  par(mar = c(4.5, 4.5, 3, 1), cex.main = 0.9)
  plot(t_seq, S_true, type = "l", lty = 1, col = "red", lwd = 3,
       xlab = "t", ylab = "S(t)",
       main = bquote("Survival function" ~ (p == .(p) * "," ~ delta[1] == .(d1) * "," ~ delta[2] == .(d2))),
       ylim = c(0, 1))
  lines(t_seq, S_bayes, lty = 2, col = "blue", lwd = 2)
  lines(t_seq, S_mle, lty = 3, col = "darkgreen", lwd = 2)
  legend("topright", c("True", "Bayes", "MLE"),
         col = c("red", "blue", "darkgreen"), lty = c(1, 2, 3), lwd = 2, bty = "n")
}

show_and_save("figure1_estimates", 900, 330 * length(Ks), fig1)
show_and_save("figure2_running_mean", 900, 330 * length(Ks), fig2)
show_and_save("figure3_ET_vs_p", 900, 600, fig3)
show_and_save("figure_survival", 900, 600, fig_surv)

# =============================================================================
# Save results
# =============================================================================
cat("\n\nTable 1:\n"); print(round(tab1, 5), row.names = FALSE)
cat("\nTable 2:\n"); print(round(tab2, 3), row.names = FALSE)
cat("\nTable 3:\n"); print(round(tab3, 5))
cat("\nTable 4:\n"); print(round(tab4_paper, 4), row.names = FALSE)

write.csv(tab1, file.path(OUT_DIR, "table1.csv"), row.names = FALSE)
write.csv(tab2, file.path(OUT_DIR, "table2.csv"), row.names = FALSE)
write.csv(tab3, file.path(OUT_DIR, "table3.csv"))
write.csv(tab4, file.path(OUT_DIR, "table4_survival_full.csv"), row.names = FALSE)
write.csv(tab4_paper, file.path(OUT_DIR, "table4_survival_paper.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))