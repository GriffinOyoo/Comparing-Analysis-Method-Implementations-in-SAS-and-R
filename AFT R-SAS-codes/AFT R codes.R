# This file contain R codes for six AFT models and their model diagnostics plots



# 1. LOG-LOGISTIC MODEL*****************************************
library(flexsurv)
(.LL <- flexsurvreg(Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR +  BMI + AFB ,
                    data = whas, dist="llogis"))
BIC(.LL)


library("eha") #When param= default"lifeAcc"
af<-aftreg(Surv(LENFOL, FSTAT)~AGE + GENDER + HR +  BMI + AFB , data = whas, dist = "loglogistic")
af
AIC(af)
BIC(af)

library(survival)
logl <- survreg(
  Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "loglogistic")
summary(logl)
AIC(logl)
BIC(logl)




# 2. LOG-NORMAL MODEL*******************************************
library(flexsurv)
(.LN <- flexsurvreg(Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,
                    data = whas, dist="lognormal"))
BIC(.LN)

library("eha")#When param= default"lifeAcc"
(af<-aftreg(Surv(LENFOL, FSTAT)~AGE + GENDER + HR +  BMI + AFB , data = whas, dist = "lognormal"))
AIC(af)
BIC(af)

library(survival)
logn <- survreg(
  Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "lognormal")
summary(logn)
AIC(logn)
BIC(logn)




# 3. WEIBULL MODEL**********************************************
library(flexsurv)
Wb <- flexsurvreg(Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,
                  data = whas, dist="weibull")
Wb
BIC(Wb)

library("eha") #When param= default"lifeAcc"
(ln<-aftreg(Surv(LENFOL, FSTAT)~AGE + GENDER + HR +  BMI + AFB , data = whas, dist = "weibull"))
AIC(ln)
BIC(ln)

library(survival)
weib <- survreg(Surv(LENFOL, FSTAT) ~AGE + GENDER + HR +  BMI + AFB,data = whas,
                dist = "weibull")
summary(weib)
AIC(weib)
BIC(weib)




# 4. EXPONENTIAL MODEL******************************************
library(flexsurv)
Wb <- flexsurvreg(Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,
                  data = whas, dist="exponential")
Wb
BIC(Wb)

library("eha")
(ln<-aftreg(Surv(LENFOL, FSTAT)~AGE + GENDER + HR +  BMI + AFB ,shape = 1, data = whas, dist = "weibull"))
AIC(ln)
BIC(ln)

library(survival)
weib <- survreg(Surv(LENFOL, FSTAT) ~AGE + GENDER + HR +  BMI + AFB,data = whas,
                dist = "exponential")
summary(weib)
AIC(weib)
BIC(weib)




# 5. GENERALIZED GAMMA MODEL************************************
library(flexsurv)
(.GG <- flexsurvreg(Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB ,
                    data = whas, dist="gengamma"))
BIC(.GG)




# 6. GENERALIZED F MODEL****************************************

library(flexsurv)
(.GF <- flexsurvreg(Surv(LENFOL, FSTAT) ~AGE + GENDER + HR + BMI + AFB,
                    data = whas, dist="genf"))
BIC(.GF)







# MODEL DIAGNOSTICS*******************************************

# 1. Standardized residuals
library(rms)
#par(mfrow = c(2, 2))
par(mfrow = c(2, 2), mar = c(4, 4, 5, 1))
gof.Weib <- psm(Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI +AFB,data = whas, dist="weibull", y=TRUE)
res.Weib <- resid(gof.Weib, type="cens")
survplot(npsurv(res.Weib ~ 1),
         conf="none",
         col="red",
         ylab="Survival probability (Gumbel)",
         xlab="Residual")
mtext("Weibull", side = 3, line = 1)
# title("Weibull")
lines(res.Weib, col="blue", lwd=0.1)

res.logn <- resid(gof.logn, type="cens")
survplot(npsurv(res.logn ~ 1),
         conf="none",
         col="red",
         main="Cox-Snell Residual Plot",
         ylab="Survival probability (Normal)",
         xlab="Residual")
mtext("Log-Normal", side = 3, line = 1)
#title("Log-Normal")
lines(res.logn, col="blue", lwd=0.1)

res.logl <- resid(gof.logl, type="cens")
survplot(npsurv(res.logl ~ 1),
         conf="none",
         col="red",
         main="Cox-Snell Residual Plot",
         ylab="Survival probability (Log-logistic)",
         xlab="Residual")
mtext("Log-Logistic", side = 3, line = 1)
#title=("Log-logistic")
lines(res.logl, col="blue", lwd=0.1)

res.exp <- resid(gof.exp, type="cens")
survplot(npsurv(res.exp ~ 1),
         conf="none",
         col="red",
         main="Cox-Snell Residual Plot",
         ylab="Survival probability (Exponential)",
         xlab="Residual")
mtext("Exponential", side = 3, line = 1)
# title="Exponential"
lines(res.exp, col="blue", lwd=0.1)





# 2.Cox-Snell residuals
library(survival)
library(flexsurv)
par(mfrow = c(2, 3)) 
fitg <- flexsurvreg(
  formula = Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "weibull")
cs <- coxsnell_flexsurvreg(fitg)
surv <- survfit(Surv(cs$est, whas$FSTAT) ~ 1)
plot(surv, fun = "cumhaz",
     xlab = "Cox–Snell residual",
     ylab = "Cumulative hazard of residuals",
     main = "Weibull")
abline(0, 1, col = "red", lwd = 2)


fitg <- flexsurvreg(
  formula = Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "llogis")
cs <- coxsnell_flexsurvreg(fitg)
surv <- survfit(Surv(cs$est, whas$FSTAT) ~ 1)
plot(surv, fun = "cumhaz",
     xlab = "Cox–Snell residual",
     ylab = "Cumulative hazard of residuals",
     main="log-logistic")
abline(0, 1, col = "red", lwd = 2)

fitg <- flexsurvreg(
  formula = Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "lognormal")
cs <- coxsnell_flexsurvreg(fitg)
surv <- survfit(Surv(cs$est, whas$FSTAT) ~ 1)
plot(surv, fun = "cumhaz",
     xlab = "Cox–Snell residual",
     ylab = "Cumulative hazard of residuals",
     main="log-normal")
abline(0, 1, col = "red", lwd = 2)

fitg <- flexsurvreg(
  formula = Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "genf")
cs <- coxsnell_flexsurvreg(fitg)
surv <- survfit(Surv(cs$est, whas$FSTAT) ~ 1)
plot(surv, fun = "cumhaz",
     xlab = "Cox–Snell residual",
     ylab = "Cumulative hazard of residuals",
     main="Generalized F")
abline(0, 1, col = "red", lwd = 2)

fitg <- flexsurvreg(
  formula = Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "gengamma")
cs <- coxsnell_flexsurvreg(fitg)
surv <- survfit(Surv(cs$est, whas$FSTAT) ~ 1)
plot(surv, fun = "cumhaz",
     xlab = "Cox–Snell residual",
     ylab = "Cumulative hazard of residuals",
     main="Generalized Gamma")
abline(0, 1, col = "red", lwd = 2)

fitg <- flexsurvreg(
  formula = Surv(LENFOL, FSTAT) ~ AGE + GENDER + HR + BMI + AFB,data = whas,
  dist = "exponential")
cs <- coxsnell_flexsurvreg(fitg)
surv <- survfit(Surv(cs$est, whas$FSTAT) ~ 1)
plot(surv, fun = "cumhaz",
     xlab = "Cox–Snell residual",
     ylab = "Cumulative hazard of residuals",
     main="Exponential")
abline(0, 1, col = "red", lwd = 2)
