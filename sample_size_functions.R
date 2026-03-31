# ==============================================================================
# Sample size calculation functions for clinical trial designs
# Endpoints: continuous, binary, time-to-event
# Designs:   superiority, non-inferiority, equivalence (TOST)
# ==============================================================================


# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------

# Critical z-value for the alpha.
# alpha : significance level (e.g. 0.025)
# sided : 1 for one-sided test, 2 for two-sided test
za <- function(alpha, sided = 2) qnorm(1 - alpha / sided)

# Critical z-value for the power (beta).
# power : desired power (e.g. 0.80)
zb <- function(power) qnorm(power)

# Pooled variance of the risk difference for binary endpoints.
# pT : event proportion in treatment group
# pC : event proportion in control group
# k  : allocation ratio (nT / nC)
.var_prop <- function(pT, pC, k) pT * (1 - pT) / k + pC * (1 - pC)

# Event probability under an exponential survival model with uniform accrual.
# Returns the probability that a patient experiences the event during the trial.
# lambda : constant hazard rate for the group
# Ta     : accrual period (same time units as lambda)
# Tb     : follow-up period after end of accrual (same time units as lambda)
prob_event_exp <- function(lambda, Ta, Tb) {
  1 - (exp(-lambda * Tb) - exp(-lambda * (Ta + Tb))) / (lambda * Ta)
}

# Core Schoenfeld-based sample size calcilation for time-to-event endpoints.
# Converts a required event count D into per-arm sample sizes accounting
# for the average event probability and dropout.
# z_num      : combined z statistic (z_alpha + z_beta)
# log_hr_adj : effective log-HR signal (varies by design; see tte_* functions)
# k          : allocation ratio (nT / nC)
# pC         : event probability in control group
# pT         : event probability in treatment group
# drop       : dropout proportion (0 = no dropout)
# Returns    : list(n_control, n_treatment, n_total, events_required)

.tte_n <- function(z_num, log_hr_adj, k, pC, pT, drop) {
  pi_avg <- (pC + pT) / 2
  D  <- ((k + 1) / k) * (z_num / log_hr_adj)^2
  nC <- D / ((k + 1) * pi_avg)
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT,
       events_required = ceiling(D))
}


# ------------------------------------------------------------------------------
# SUPERIORITY
# ------------------------------------------------------------------------------

# Sample size for a superiority trial with a continuous endpoint.
# alpha : one-sided significance level (default 0.025)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1 = equal allocation)
# SD    : pooled standard deviation of the outcome
# d     : expected mean difference (mu_T - mu_C)
# delta : superiority margin (must be < d)
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total)

cont_super <- function(alpha = 0.025, power = 0.80, k = 1,
                       SD, d, delta, drop = 0) {
  z  <- za(alpha, sided = 1) + zb(power)
  nC <- (1 + 1/k) * SD^2 * (z / (d - delta))^2
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT)
}

# Sample size for a superiority trial with a binary endpoint.
# alpha : one-sided significance level (default 0.025)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# pT    : expected event proportion in treatment group
# pC    : expected event proportion in control group
# delta : superiority margin on the risk-difference scale (must be < pT - pC)
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total)

prop_super <- function(alpha = 0.025, power = 0.80, k = 1,
                       pT, pC, delta, drop = 0) {
  d  <- pT - pC
  z  <- za(alpha, sided = 1) + zb(power)
  nC <- (z / (d - delta))^2 * .var_prop(pT, pC, k)
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT)
}

# Sample size for a superiority trial with a time-to-event endpoint.
# alpha : one-sided significance level (default 0.025)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# HR    : expected hazard ratio (treatment / control); HR > 1 means faster event
# pC    : event probability in control group (use prob_event_exp)
# pT    : event probability in treatment group (use prob_event_exp)
# delta : superiority margin on the log-HR scale (0 = pure superiority)
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total, events_required)
tte_super <- function(alpha = 0.025, power = 0.80, k = 1,
                      HR, pC, pT, delta, drop = 0) {
  z   <- za(alpha, sided = 1) + zb(power)
  adj <- log(HR) - delta
  .tte_n(z, adj, k, pC, pT, drop)
}


# ------------------------------------------------------------------------------
# NON-INFERIORITY
# ------------------------------------------------------------------------------

# Sample size for a non-inferiority trial with a continuous endpoint.
# alpha : one-sided significance level (default 0.025)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# SD    : pooled standard deviation of the outcome
# d     : assumed true mean difference (default 0 = drugs perform equally)
# delta : non-inferiority margin (maximum acceptable inferiority)
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total)

cont_noninf <- function(alpha = 0.025, power = 0.80, k = 1,
                        SD, d = 0, delta, drop = 0) {
  z  <- za(alpha, sided = 1) + zb(power)
  nC <- (1 + 1/k) * SD^2 * (z / (d + delta))^2
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT)
}

# Sample size for a non-inferiority trial with a binary endpoint.

# alpha : one-sided significance level (default 0.025)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# pT    : expected event proportion in treatment group
# pC    : expected event proportion in control group
# delta : non-inferiority margin on the risk-difference scale
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total)

prop_noninf <- function(alpha = 0.025, power = 0.80, k = 1,
                        pT, pC, delta, drop = 0) {
  d  <- pT - pC
  z  <- za(alpha, sided = 1) + zb(power)
  nC <- (z / (d + delta))^2 * .var_prop(pT, pC, k)
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT)
}

# Sample size for a non-inferiority trial with a time-to-event endpoint.
# alpha : one-sided significance level (default 0.025)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# HR    : assumed true hazard ratio (1.0 = drugs perform equally)
# pC    : event probability in control group (use prob_event_exp)
# pT    : event probability in treatment group (use prob_event_exp)
# delta : non-inferiority margin on the log-HR scale (e.g. log(1.3))
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total, events_required)

tte_noninf <- function(alpha = 0.025, power = 0.80, k = 1,
                       HR, pC, pT, delta, drop = 0) {
  z   <- za(alpha, sided = 1) + zb(power)
  adj <- log(HR) + delta
  .tte_n(z, adj, k, pC, pT, drop)
}


# ------------------------------------------------------------------------------
# EQUIVALENCE (TOST)
# ------------------------------------------------------------------------------

# Sample size for an equivalence trial with a continuous endpoint (TOST).
# alpha : one-sided significance level per TOST component (default 0.05)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# SD    : pooled standard deviation of the outcome
# d     : assumed true mean difference (default 0 = drugs identical)
# delta : equivalence margin; drugs are equivalent if |d| < delta
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total)

cont_equiv <- function(alpha = 0.05, power = 0.80, k = 1,
                       SD, d = 0, delta, drop = 0) {
  z  <- za(alpha, sided = 1) + zb(power)
  nC <- (1 + 1/k) * SD^2 * (z / (delta - abs(d)))^2
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT)
}

# Sample size for an equivalence trial with a binary endpoint (TOST).
# alpha : one-sided significance level per TOST component (default 0.05)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# pT    : expected event proportion in treatment group
# pC    : expected event proportion in control group
# delta : equivalence margin on the risk-difference scale
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total)

prop_equiv <- function(alpha = 0.05, power = 0.80, k = 1,
                       pT, pC, delta, drop = 0) {
  d  <- pT - pC
  z  <- za(alpha, sided = 1) + zb(power)
  nC <- (z / (delta - abs(d)))^2 * .var_prop(pT, pC, k)
  nC <- ceiling(nC / (1 - drop))
  nT <- ceiling(k * nC)
  list(n_control = nC, n_treatment = nT, n_total = nC + nT)
}

# Sample size for an equivalence trial with a time-to-event endpoint (TOST).
# alpha : one-sided significance level per TOST component (default 0.05)
# power : desired power (default 0.80)
# k     : allocation ratio nT / nC (default 1)
# HR    : assumed true hazard ratio (1.0 = drugs identical)
# pC    : event probability in control group (use prob_event_exp)
# pT    : event probability in treatment group (use prob_event_exp)
# delta : equivalence margin on the log-HR scale (e.g. log(1.3))
# drop  : dropout proportion (default 0)
# Returns: list(n_control, n_treatment, n_total, events_required)
tte_equiv <- function(alpha = 0.05, power = 0.80, k = 1,
                      HR, pC, pT, delta, drop = 0) {
  z   <- za(alpha, sided = 1) + zb(power)
  adj <- delta - abs(log(HR))
  .tte_n(z, adj, k, pC, pT, drop)
}
