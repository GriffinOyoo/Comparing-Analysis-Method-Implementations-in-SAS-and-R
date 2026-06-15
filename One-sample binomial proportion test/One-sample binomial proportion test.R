# H0: P=0.5#############################################

# 1. Exact Binomial Test.
binom.test(heads_count, total_flips, p = 0.5, conf.level = 0.95)

# 2. Wald(Asymptotic) Test.
library(DescTools)
p0 <- 0.5
phat <- heads_count / total_flips
z <- (phat - p0) / sqrt(phat * (1 - phat) / total_flips)
p_two_sided <- 2 * (1 - pnorm(abs(z)))
p_two_sided
BinomCI(heads_count, total_flips, method = "wald")

# 3. Mid-P adjusted Exact Binomial Test.
library(exactci)
binom.exact(heads_count,total_flips, p = 0.5,alternative = "greater", midp = TRUE,tsmethod = "central")
binom.exact(heads_count,total_flips, p = 0.5,alternative = "less", midp = TRUE,tsmethod = "central")
binom.exact(heads_count,total_flips, p = 0.5, midp = TRUE,tsmethod = "central")

# 3. Wilson Score Test.
prop.test(heads_count, total_flips, p = 0.5, correct = FALSE)



#H_0 : p = 0.19########################################

# 1. Exact Binomial Test.
binom.test(heads_count, total_flips, p = 0.19, conf.level = 0.95)

# 2. Wald(Asymptotic) Test.
library(DescTools)
p0 <- 0.19
phat <- heads_count / total_flips
z <- (phat - p0) / sqrt(phat * (1 - phat) / total_flips)
p_two_sided <- 2 * (1 - pnorm(abs(z)))
p_two_sided
BinomCI(heads_count, total_flips, method = "wald")

# 3. Mid-P adjusted Exact Binomial Test.
library(exactci)
binom.exact(heads_count,total_flips, p = 0.19,alternative = "greater", midp = TRUE,tsmethod = "central")
binom.exact(heads_count,total_flips, p = 0.19,alternative = "less", midp = TRUE,tsmethod = "central")
binom.exact(heads_count,total_flips, p = 0.19, midp = TRUE,tsmethod = "central")

# 3. Wilson Score Test.
prop.test(heads_count, total_flips, p = 0.19, correct = FALSE)
