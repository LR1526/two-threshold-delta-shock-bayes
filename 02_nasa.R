# =============================================================================
# NASA C-MAPSS (FD001) analysis for
# "Bayesian Inference for Survival Analysis in Two-Threshold Generalized
#  delta-shock Models"  -- revision, Reviewer 1, comment 7.
#
# Produces: Table 4 (estimates of theta), Table 5 (survival function),
#           Table 6 (goodness of fit: Pearson + split predictive check),
#           Table 7 (sensitivity of the fit to the choice of mapping),
#           Table 8 (repeated 20-engine subsamples).
#
# R >= 4.0, base R only. Put train_FD001.txt in the working directory.
# =============================================================================

set.seed(2026)
OUT_DIR   <- "outputs_nasa"
DATA_FILE <- "train_FD001.txt"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------- settings --
SENSOR <- "s12"      # degradation indicator (sensor 12)
Z1     <- 1.0        # a cycle is a SHOCK if |standardised residual| > Z1
Z2     <- 1.5        # a shock is SEVERE  if |standardised residual| > Z2
d1     <- 1          # deterministic region: inter-shock time <= d1
d2     <- 3          # probabilistic region: d1 < inter-shock time <= d2
N_SMALL   <- 20      # small sample (engines 1..20, deterministic)
R_SUB     <- 200     # random 20-engine subsamples (Table 8)
N_SPLITS  <- 200     # random splits for the split predictive check
NREP_SPC  <- 100     # posterior draws per split
RUN_GRID  <- TRUE    # Table 7 (sensitivity to mapping)
TS        <- c(5, 10, 20)

# ------------------------------------------------------------------- data ---
dat <- read.table(DATA_FILE, header = FALSE)
names(dat) <- c("unit", "cycle", "os1", "os2", "os3", paste0("s", 1:21))
units <- sort(unique(dat$unit))

# per-engine standardised residuals after removing a cubic trend (the slow
# degradation drift); what remains is the cycle-to-cycle fluctuation.
std_resid <- function(x) {
  tt <- seq_along(x) / length(x)
  r  <- resid(lm(x ~ poly(tt, 3)))
  r / sqrt(mean(r^2))
}
get_Z <- function(sensor) {
  lapply(units, function(u) std_resid(dat[dat$unit == u, sensor]))
}

# ---------------------------------------------------------------------------
# MAPPING (inter-shock times are the X_i of the model, NOT sensor levels)
#   shock      : cycle t with |z_t| > z1
#   X_i        : number of cycles between consecutive shocks (integer, >= 1)
#   region     : X <= d1 deterministic; d1 < X <= d2 probabilistic; X > d2 harmless
#   probabilistic region: the shock is FATAL iff it is severe (|z| > z2);
#                         theta = P(severe | shock arriving after d1 < X <= d2)
#   lifetime T : cycle of the first fatal shock (one lifetime per engine);
#                lifetimes only are used for inference, as in Section 3.1.
# ---------------------------------------------------------------------------
first_T <- function(z, z1, z2) {
  idx <- which(abs(z) > z1)
  if (length(idx) < 2) return(NA_real_)
  gaps  <- diff(idx)
  later <- idx[-1]
  fail  <- gaps <= d1 | (gaps <= d2 & abs(z[later]) > z2)
  if (!any(fail)) return(NA_real_)
  later[which(fail)[1]]
}

# all (renewal) lifetimes of an engine: after a fatal shock the process restarts
episode_T <- function(z, z1, z2) {
  idx <- which(abs(z) > z1)
  out <- numeric(0); r <- 0; prev <- NA
  for (s in idx) {
    if (is.na(prev)) { prev <- s; next }
    x <- s - prev
    if (x <= d1 || (x <= d2 && abs(z[s]) > z2)) {
      out <- c(out, s - r); r <- s; prev <- NA
    } else prev <- s
  }
  out
}

prepare_mapping <- function(sensor, z1, z2) {
  Z   <- get_Z(sensor)
  T1  <- sapply(Z, first_T, z1 = z1, z2 = z2)
  ep  <- unlist(lapply(Z, episode_T, z1 = z1, z2 = z2))
  nsh <- sum(sapply(Z, function(z) sum(abs(z) > z1)))
  ncy <- sum(sapply(Z, length))
  mid_n <- 0; mid_s <- 0
  for (z in Z) {
    idx <- which(abs(z) > z1)
    if (length(idx) > 1) {
      g <- diff(idx); later <- idx[-1]
      m <- g > d1 & g <= d2
      mid_n <- mid_n + sum(m)
      mid_s <- mid_s + sum(abs(z[later[m]]) > z2)
    }
  }
  list(T1 = T1, ep = ep, p_hat = round(nsh / ncy, 3),
       theta_ref = mid_s / mid_n, n_mid = mid_n, n_cens = sum(is.na(T1)), Z = Z)
}

# ------------------------------------------------- model (same as paper code)
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

set_model <- function(p_, tmax_) {
  p      <<- p_
  TMAX   <<- tmax_
  CC     <<- setup_C(p_, d1, d2, tmax_)
  Dcache <<- vector("list", tmax_)
}

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
  g <- seq(0.001, 0.999, length.out = 2000)
  l <- loglik(P, g)
  lmax <- loglik(P, th_hat)
  ok <- g[2 * (lmax - l) <= qchisq(0.95, 1)]
  if (length(ok) == 0) return(c(0.001, 0.999))
  c(min(ok), max(ok))
}

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

pmf_theta <- function(x) {
  m <- 0:(ncol(CC$C1) - 1)
  as.vector((x * CC$C1 + CC$C2) %*% (1 - x)^m)
}

sample_theta <- function(P, k, al = 1, be = 1) {
  w <- post_w(P, al, be)
  i <- sample.int(length(w), k, TRUE, w)
  rbeta(k, al + P$a[i], be + P$m[i])
}

analyse <- function(tt) {
  P  <- make_post(coef_array(tt))
  th <- mle(P)
  list(P = P, mle = th, bayes = post_mean(P), cr = post_ci(P), lr = lr_ci(P, th))
}

# ------------------------------------------------ goodness-of-fit tools -----
pearson_gof <- function(tt, th) {
  f <- pmf_theta(th)[1:TMAX]
  obs <- tabulate(tt, TMAX)
  E <- length(tt) * f
  ob <- numeric(0); ex <- numeric(0); co <- 0; ce <- 0
  for (i in seq_len(TMAX)) {
    co <- co + obs[i]; ce <- ce + E[i]
    if (ce >= 5) { ob <- c(ob, co); ex <- c(ex, ce); co <- 0; ce <- 0 }
  }
  if (ce > 0) { k <- length(ob); ob[k] <- ob[k] + co; ex[k] <- ex[k] + ce }
  X2 <- sum((ob - ex)^2 / ex)
  dfree <- length(ob) - 2
  c(X2 = X2, df = dfree, p = pchisq(X2, dfree, lower.tail = FALSE))
}

spc_p <- function(tt, stat = sd, q = 0.5, nrep = NREP_SPC) {
  N <- length(tt)
  nobs <- ceiling(q * N)
  idx <- sample.int(N)
  obs <- tt[idx[1:nobs]]
  ho  <- tt[idx[(nobs + 1):N]]
  th  <- sample_theta(make_post(coef_array(obs)), nrep)
  d_new <- stat(ho)
  cnt <- 0
  for (x in th) {
    yrep <- sample.int(TMAX, length(ho), TRUE, pmf_theta(x))
    d_rep <- stat(yrep)
    cnt <- cnt + (d_rep > d_new) + 0.5 * (d_rep == d_new)
  }
  pp <- cnt / nrep
  2 * min(pp, 1 - pp)
}
q90 <- function(x) quantile(x, 0.9, names = FALSE)

# =============================================================================
# MAIN MAPPING
# =============================================================================
mp <- prepare_mapping(SENSOR, Z1, Z2)
stopifnot(mp$n_cens == 0)
T100 <- mp$T1
TM <- max(c(mp$T1, mp$ep)) + 5
set_model(mp$p_hat, TM)
th_ref <- mp$theta_ref
S_ref  <- function(t) mean(mp$ep > t)

cat(sprintf("\nMapping: sensor %s, shock |z|>%.1f, severe |z|>%.1f, delta1=%d, delta2=%d\n",
            SENSOR, Z1, Z2, d1, d2))
cat(sprintf("p_hat = %.3f | n = %d first-failure lifetimes (mean %.1f, sd %.1f)\n",
            mp$p_hat, length(T100), mean(T100), sd(T100)))
cat(sprintf("reference theta (all %d intermediate-gap shocks) = %.4f\n", mp$n_mid, th_ref))
cat(sprintf("reference survival from %d pooled episodes\n", length(mp$ep)))

for (th in c(0.1, 0.5, 0.9)) {
  m <- 0:(ncol(CC$C1) - 1)
  f <- sapply(1:TMAX, function(t) sum((1 - th)^m * (th * CC$C1[t, ] + CC$C2[t, ])))
  stopifnot(sum(f) <= 1 + 1e-8, sum(f) > 0.99)
}

# ---------------------------------------------------------- Tables 4 and 5 --
a100 <- analyse(T100)
a20  <- analyse(T100[1:N_SMALL])
A <- list("Large sample (100 engines)" = a100,
          "Small sample (20 engines)"  = a20)

tab4 <- data.frame()
for (sc in names(A)) {
  a <- A[[sc]]
  tab4 <- rbind(tab4,
    data.frame(Scenario = sc, Method = "Bayesian (Beta(1,1))", theta_hat = a$bayes,
               lower = a$cr[1], upper = a$cr[2], width = diff(a$cr), theta_ref = th_ref),
    data.frame(Scenario = sc, Method = "MLE (LR interval)", theta_hat = a$mle,
               lower = a$lr[1], upper = a$lr[2], width = diff(a$lr), theta_ref = th_ref))
}

tab5 <- data.frame()
for (t0 in TS) {
  Sf <- make_S(t0)
  for (sc in names(A)) {
    a <- A[[sc]]
    sb <- post_S_mean(a$P, t0)
    sm <- Sf(a$mle)
    tab5 <- rbind(tab5,
      data.frame(t = t0, Scenario = sc, Method = "Bayesian", S_hat = sb,
                 S_ref = S_ref(t0), abs_error = abs(sb - S_ref(t0))),
      data.frame(t = t0, Scenario = sc, Method = "MLE", S_hat = sm,
                 S_ref = S_ref(t0), abs_error = abs(sm - S_ref(t0))))
  }
}

cat("\nTable 4: estimation of theta\n"); print(format(tab4, digits = 3), row.names = FALSE)
cat("\nTable 5: survival function S(t)\n"); print(format(tab5, digits = 3), row.names = FALSE)

# ------------------------------------------------------ Table 6: model fit --
cat("\nTable 6: goodness of fit (100 engines)\n")
pg <- pearson_gof(T100, a100$mle)
ps_sd  <- replicate(N_SPLITS, spc_p(T100, sd))
ps_q90 <- replicate(N_SPLITS, spc_p(T100, q90))
tab6 <- data.frame(
  check = c("Pearson chi-square", "Split predictive check, T = sd",
            "Split predictive check, T = 90th percentile"),
  statistic = c(pg["X2"], NA, NA),
  df = c(pg["df"], NA, NA),
  p_value_or_median_p = c(pg["p"], median(ps_sd), median(ps_q90)),
  rejection_rate_at_0.05 = c(NA, mean(ps_sd < 0.05), mean(ps_q90 < 0.05)))
print(format(tab6, digits = 3), row.names = FALSE)

# negative control
Lfull <- as.numeric(tapply(dat$cycle, dat$unit, max))
meanT_p <- function(pp, th) {
  qq <- 1 - pp; pp2 <- qq^d1 - qq^d2; pp3 <- 1 - qq^d1; u <- pp2 * th + pp3
  (1 + u) / (pp * u)
}
sim_life <- function(theta, pp) {
  t <- rgeom(1, pp) + 1
  repeat {
    x <- rgeom(1, pp) + 1
    t <- t + x
    if (x <= d1) break
    if (x <= d2 && runif(1) < theta) break
  }
  t
}
neg <- data.frame()
for (th in c(0.1, 0.5, 0.9)) {
  pm <- uniroot(function(pp) meanT_p(pp, th) - mean(Lfull), c(1e-4, 0.9))$root
  sds <- replicate(300, sd(replicate(length(Lfull), sim_life(th, pm))))
  neg <- rbind(neg, data.frame(theta = th, p_matched = pm, obs_sd = sd(Lfull),
                               model_sd_median = median(sds),
                               P_model_sd_le_obs = mean(sds <= sd(Lfull))))
}
cat("\nNegative control (whole-engine lifetimes):\n"); print(format(neg, digits = 3), row.names = FALSE)

# ------------------------- Table 8: repeated random 20-engine subsamples ----
res <- matrix(NA, R_SUB, 8)
for (r in 1:R_SUB) {
  a <- analyse(T100[sample.int(length(T100), N_SMALL)])
  res[r, ] <- c(a$bayes, a$mle,
                a$cr[1] <= th_ref & th_ref <= a$cr[2],
                a$lr[1] <= th_ref & th_ref <= a$lr[2],
                diff(a$cr), diff(a$lr),
                post_S_mean(a$P, 10), make_S(10)(a$mle))
}
tab8 <- data.frame(
  method = c("Bayesian", "MLE"),
  bias = c(mean(res[, 1]), mean(res[, 2])) - th_ref,
  MSE = c(mean((res[, 1] - th_ref)^2), mean((res[, 2] - th_ref)^2)),
  coverage = c(mean(res[, 3]), mean(res[, 4])),
  mean_width = c(mean(res[, 5]), mean(res[, 6])),
  MSE_S10 = c(mean((res[, 7] - S_ref(10))^2), mean((res[, 8] - S_ref(10))^2)))
cat("\nTable 8: repeated 20-engine subsamples (R =", R_SUB, ")\n")
print(format(tab8, digits = 3), row.names = FALSE)

# ---------------------------- Table 7: sensitivity to the choice of mapping --
tab7 <- data.frame()
if (RUN_GRID) {
  for (sn in c("s12", "s11", "s7", "s20", "s8", "s13")) {
    for (zz in list(c(1.0, 1.5), c(1.2, 1.8))) {
      m2 <- prepare_mapping(sn, zz[1], zz[2])
      if (m2$n_cens > 0) next
      set_model(m2$p_hat, max(c(m2$T1, m2$ep)) + 5)
      a <- analyse(m2$T1)
      pr <- replicate(50, spc_p(m2$T1, sd))
      tab7 <- rbind(tab7, data.frame(
        sensor = sn, z1 = zz[1], z2 = zz[2], p_hat = m2$p_hat,
        theta_ref = m2$theta_ref, theta_MLE = a$mle, theta_Bayes = a$bayes,
        SPC_median_p = median(pr), SPC_reject_rate = mean(pr < 0.05)))
    }
  }
  cat("\nTable 7: fit under alternative mappings\n"); print(format(tab7, digits = 3), row.names = FALSE)
  set_model(mp$p_hat, TM)
}

write.csv(tab4, file.path(OUT_DIR, "table4_theta.csv"), row.names = FALSE)
write.csv(tab5, file.path(OUT_DIR, "table5_survival.csv"), row.names = FALSE)
write.csv(tab6, file.path(OUT_DIR, "table6_gof.csv"), row.names = FALSE)
write.csv(neg,  file.path(OUT_DIR, "table6b_negative_control.csv"), row.names = FALSE)
write.csv(tab7, file.path(OUT_DIR, "table7_mapping_sensitivity.csv"), row.names = FALSE)
write.csv(tab8, file.path(OUT_DIR, "table8_subsamples.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))