

#  1. WEIGHTED LOG-RANK TEST

#  PETO-PETO*********************************************************

# a) survminer::surv_pvalue
fit <- surv_fit(Surv(LENFOL, FSTAT)~AFB, data = WHAS)
surv_pvalue(fit,
            method = "S1")
# Extending the dp for the p-value
p<-surv_pvalue(fit, method="S1")$pval
sprintf("%.9f", p)



#  b) Custom R function for Peto-Peto test
peto_peto_test <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  surv <- mf[[1]]
  time <- surv[, "time"]
  status <- ifelse(surv[, "status"] == 2, 1, surv[, "status"])
  group <- as.factor(mf[[2]])
  groups <- levels(group)
  event_times <- sort(unique(time[status == 1]))
  k <- length(event_times)   
  g <- length(groups)       
  n_jk <- d_jk <- matrix(0, k, g)
  for (j in seq_len(g)) {
    t_j <- time[group == groups[j]]
    s_j <- status[group == groups[j]]
    for (i in seq_len(k)) {
      ti <- event_times[i]
      n_jk[i, j] <- sum(t_j >= ti)
      d_jk[i, j] <- sum(t_j == ti & s_j == 1)
    }
  }
  n_i <- rowSums(n_jk)
  d_i <- rowSums(d_jk)
  w <- cumprod(1 - d_i / (n_i + 1)) 
  e1 <- n_jk[, 1] * d_i / n_i
  Q <- sum(w * (d_jk[, 1] - e1))
  VarQ <- sum(
    w^2 * n_jk[, 1] * n_jk[, 2] * d_i * (n_i - d_i) /
      (n_i^2 * (n_i - 1)),
    na.rm = TRUE
  )
  chisq <- (Q^2) / VarQ
  pvalue <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  list(
    chisq = chisq,
    df = 1,
    pvalue = pvalue,
    method = "Peto-Peto (S1)"
  )
}
peto_peto_test(Surv(LENFOL, FSTAT)~AFB, data = WHAS)

# c)coin::logrank_test
#Asymptotic
WHAS$AFB <- as.factor(WHAS$AFB)
logrank_test(Surv(WHAS$LENFOL, WHAS$FSTAT) ~ AFB, data = WHAS,
             type = "Peto-Peto",
             distribution = "asymptotic")

# d) survival::survdiff
(L<- survival::survdiff(
  survival::Surv(LENFOL, FSTAT)~AFB,rho = 1,data = WHAS))

#  MODIFIED PETO-PETO*********************************************************

#  a)survminer::surv_pvalue
fit.null <- surv_fit(Surv(LENFOL,FSTAT)~AFB, data = WHAS)
surv_pvalue(fit.null,
            method = "S2")
# Extending the dp for the p-value
p<-surv_pvalue(fit, method="S2")$pval
sprintf("%.9f", p)
 

#  b) Custom R function for modified Peto-Peto test
modified_peto_peto_test <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  surv <- mf[[1]]
  time <- surv[, "time"]
  status <- ifelse(surv[, "status"] == 2, 1, surv[, "status"])
  group <- as.factor(mf[[2]])
  groups <- levels(group)
  event_times <- sort(unique(time[status == 1]))
  k <- length(event_times)   
  g <- length(groups)       
  n_jk <- d_jk <- matrix(0, k, g)
  for (j in seq_len(g)) {
    t_j <- time[group == groups[j]]
    s_j <- status[group == groups[j]]
    for (i in seq_len(k)) {
      ti <- event_times[i]
      n_jk[i, j] <- sum(t_j >= ti)
      d_jk[i, j] <- sum(t_j == ti & s_j == 1)
    }
  }
  n_i <- rowSums(n_jk)
  d_i <- rowSums(d_jk)
  s1 <- cumprod(1 - d_i / (n_i + 1))
  w  <- s1 * n_i / (n_i + 1)
  e1 <- n_jk[, 1] * d_i / n_i
  Q <- sum(w * (d_jk[, 1] - e1))
  VarQ <- sum(
    w^2 * n_jk[, 1] * n_jk[, 2] * d_i * (n_i - d_i) /
      (n_i^2 * (n_i - 1)),
    na.rm = TRUE
  )
  chisq <- (Q^2) / VarQ
  pvalue <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  list(chisq = chisq,df = 1,pvalue = pvalue,method = "Modified Peto-Peto (S2)"
  )
}
modified_peto_peto_test(Surv(LENFOL, FSTAT)~AFB, data = WHAS)

#  TARONE-WARE *********************************************************

# a) survminer::surv_pvalue
fit<-survfit(Surv(LENFOL, FSTAT)~AFB, data = WHAS)
surv_pvalue(fit,method="sqrtN")
p<-surv_pvalue(fit, method="SqrtN")$pval
sprintf("%.9f", p)

# b) Custom R function for Tarone-Ware test
tarone_ware_test <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  surv <- mf[[1]]
  time <- surv[, "time"]
  status <- ifelse(surv[, "status"] == 2, 1, surv[, "status"])
  group <- as.factor(mf[[2]])
  groups <- levels(group)
  event_times <- sort(unique(time[status == 1]))
  k <- length(event_times)
  g <- length(groups)
  n_jk <- d_jk <- matrix(0, k, g)
  for (j in seq_len(g)) {
    t_j <- time[group == groups[j]]
    s_j <- status[group == groups[j]]
    for (i in seq_len(k)) {
      ti <- event_times[i]
      n_jk[i,j] <- sum(t_j >= ti)
      d_jk[i,j] <- sum(t_j == ti & s_j == 1)
    }
  }
  n_i <- rowSums(n_jk)
  d_i <- rowSums(d_jk)
  w <- sqrt(n_i)   
  e1 <- n_jk[,1] * d_i / n_i
  Q <- sum(w * (d_jk[,1] - e1))
  VarQ <- sum(
    w^2 * n_jk[,1] * n_jk[,2] * d_i * (n_i - d_i) /
      (n_i^2 * (n_i - 1)),
    na.rm = TRUE
  )
  chisq <- (Q^2) / VarQ
  pvalue <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  list(chisq = chisq, df = 1, pvalue = pvalue, method = "Tarone-Ware")
}
tarone_ware_test(Surv(LENFOL, FSTAT)~AFB, data = WHAS)

#c) coin::logrank_test
logrank_test(Surv(LENFOL, FSTAT)~AFB,data = WHAS,type = "Tarone-Ware",
             distribution = "asymptotic"
)


#  GEHAN-BRESLOW *********************************************************

# a) survminer::surv_pvalue
fit<-survfit(Surv(LENFOL, FSTAT)~AFB, data = WHAS)
library(survminer)
surv_pvalue(fit,method="n")
p<-surv_pvalue(fit, method="n")$pval
sprintf("%.9f", p)

# b) Custom R function for Gehan-breslow test
Gehan_Breslow_test <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  surv <- mf[[1]]
  time <- surv[, "time"]
  status <- ifelse(surv[, "status"] == 2, 1, surv[, "status"])
  group <- as.factor(mf[[2]])
  groups <- levels(group)
  event_times <- sort(unique(time[status == 1]))
  k <- length(event_times)
  g <- length(groups)
  n_jk <- d_jk <- matrix(0, k, g)
  for (j in seq_len(g)) {
    t_j <- time[group == groups[j]]
    s_j <- status[group == groups[j]]
    for (i in seq_len(k)) {
      ti <- event_times[i]
      n_jk[i,j] <- sum(t_j >= ti)
      d_jk[i,j] <- sum(t_j == ti & s_j == 1)
    }
  }
  n_i <- rowSums(n_jk)
  d_i <- rowSums(d_jk)
  w <- n_i  
  e1 <- n_jk[,1] * d_i / n_i
  Q <- sum(w * (d_jk[,1] - e1))
  VarQ <- sum(
    w^2 * n_jk[,1] * n_jk[,2] * d_i * (n_i - d_i) /
      (n_i^2 * (n_i - 1)),
    na.rm = TRUE
  )
  chisq <- (Q^2) / VarQ
  pvalue <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  list(chisq = chisq, df = 1, pvalue = pvalue, method = "Gehan-Breslow")
}
Gehan_Breslow_test(Surv(LENFOL, FSTAT)~AFB, data = WHAS)

# c) coin::logrank_test
logrank_test(Surv(LENFOL, FSTAT)~AFB,data = WHAS,type = "Gehan-Breslow",
             distribution = "asymptotic"
)



#FLEMING-HARRINGTON *******************************************************

# a). nphRCT::wlrt
#Places weight on hazards at the end of the study.
(ppp <- nphRCT::wlrt(survival::Surv(LENFOL, FSTAT)~AFB,   
                     data = WHAS,method = "fh",rho = 0,gamma = 1 ))
p<-2 * (1 - pnorm(abs(ppp$z), 0, 1))
cat("chisquare:", round(ppp$z^2, 4), "  p-value:", round(p, 4), "\n")

# b) coin::logrank_test
#Places weight on hazards at the end of the study.
logrank_test(Surv(LENFOL, FSTAT)~AFB, data = WHAS,
             type = "Fleming-Harrington", rho = 0, gamma = 1,
             distribution = "asymptotic")






# 2. MAXCOMBO TEST.

# a). MaxCombo(No strata)*********************************************
set.seed(254)
logrank.maxtest(
  time   = WHAS$LENFOL,
  event  = WHAS$FSTAT,
  group  = WHAS$AFB,
  alternative = c("two.sided"),
  rho = c(0, 1, 0,1),
  gamma = c(0, 0, 1, 1),
  event_time_weights = NULL,
  algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 1e-05, releps = 0)
)
#b). Stratified MaxCombo**********************************************
SMCtest(
  time= WHAS$LENFOL,
  event= WHAS$FSTAT,
  group=WHAS$AFB,
  stratum=WHAS$GENDER,
  alternative = c("two.sided"),
  rho = c(0, 1, 0, 1),
  gamma = c(0, 0, 1, 1)
)


# 3. MODESTLY WEIGHTED LOG-RANK TEST 
# a). Modestly Weighted Log-rank test(No strata) *********************
modestly <- nphRCT::wlrt(Surv(LENFOL, FSTAT)~AFB,data = WHAS,method = "mw",
                         s_star = 0.5
                         #t_star = 450
)
chi_sq<-(modestly$z)^2
p_value<-pchisq((modestly$z)^2, df = 1, lower.tail = FALSE)
cat(sprintf("Chi-square: %.6f   p-value: %.6f\n", chi_sq, p_value))

# b). Stratified Modestly Weighted Log-Rank Test (Stratified MWLRT)******8
smodestly <- nphRCT::wlrt(
  formula = Surv(LENFOL, FSTAT)~AFB + strata(GENDER),data = WHAS,method = "mw",
  s_star = 0.5)
z_strata <- smodestly$by_strata$z
z_combined <- smodestly$combined$z
pchi_sq <- z_strata^2
cchi_sq <- z_combined^2
pp_value <- pchisq(pchi_sq, df = 1, lower.tail = FALSE)
cp_value <- pchisq(cchi_sq, df = 1, lower.tail = FALSE)
cat(sprintf("\nStrata Chi-square: %.6f   p-value: %.6f\n", pchi_sq, pp_value))
cat(sprintf("\nCombined Chi-square: %.6f   p-value: %.6f\n", cchi_sq, cp_value))
