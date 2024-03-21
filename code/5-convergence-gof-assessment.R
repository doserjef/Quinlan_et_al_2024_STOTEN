# 5-convergence-gof-assessment.R: assess model convergence and goodness of 
#                                 fit for all candidate models
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)

# Load data ---------------------------------------------------------------
load("data/spAbundance-data.rda")
sp.names <- dimnames(data.list$y)[[1]]

# Read in all models ------------------------------------------------------
# Each load will read in the following objects for each species.
#    beta.samples -> the spcies-specific regression coefficients
#    beta.star.samples -> the species-specific random day effects
#    lambda.samples -> the factor loadings, not all that releveant
#    fit.y -> test statistic values for real data from a Bayesian PPC
#    fit.y.rep -> test statistic values for fake data from a Bayesian PPC
#    kappa.samples -> NB dispersion values (lower = more overdispersion)
#    y.rep.means -> mean fitted value for each species/site/rep
#    waic.model -> results from WAIC
#    log.like.by.sp -> log likelihood values by species
# Null only ---------------------------
load("results/lfMsAbund-null-samples.rda")
# GoF assessment
# Plot of true vs. fitted for each species. Looks pretty good for a NB model
for (i in 1:length(sp.names)) {
  plot(c(data.list$y[i, , ]), c(y.rep.means[i, , ]), pch = 19)
  abline(0, 1)
  Sys.sleep(1)
}
# Bayesian p-value for each species
N <- nrow(data.list$y)
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}
hist(bpvs, xlim = c(0, 1))
# Convergence assessment
ess.beta <- effectiveSize(beta.samples)
summary(ess.beta)
plot(beta.samples, density = FALSE)
ess.beta.star <- effectiveSize(beta.star.samples)
summary(ess.beta.star)
plot(kappa.samples, density = FALSE)
ess.kappa <- effectiveSize(kappa.samples)
ess.lambda <- effectiveSize(lambda.samples)
summary(ess.lambda[ess.lambda != 0])
plot(lambda.samples, density = FALSE)
# Devel only --------------------------
load("results/lfMsAbund-devel-samples.rda")
# GoF assessment
# Plot of true vs. fitted for each species. Looks pretty good for a NB model
for (i in 1:length(sp.names)) {
  plot(c(data.list$y[i, , ]), c(y.rep.means[i, , ]), pch = 19)
  abline(0, 1)
  Sys.sleep(1)
}
# Bayesian p-value for each species
N <- nrow(data.list$y)
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}
hist(bpvs, xlim = c(0, 1))
# Convergence assessment
ess.beta <- effectiveSize(beta.samples)
summary(ess.beta)
# plot(beta.samples, density = FALSE)
ess.beta.star <- effectiveSize(beta.star.samples)
summary(ess.beta.star)
# plot(kappa.samples, density = FALSE)
ess.kappa <- effectiveSize(kappa.samples)
summary(ess.kappa)
ess.lambda <- effectiveSize(lambda.samples)
summary(ess.lambda[ess.lambda != 0])
# plot(lambda.samples, density = FALSE)

# Apiary only --------------------------
load("results/lfMsAbund-apiary-samples.rda")
# GoF assessment
# Plot of true vs. fitted for each species. Looks pretty good for a NB model
for (i in 1:length(sp.names)) {
  plot(c(data.list$y[i, , ]), c(y.rep.means[i, , ]), pch = 19)
  abline(0, 1)
  Sys.sleep(1)
}
# Bayesian p-value for each species
N <- nrow(data.list$y)
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}
hist(bpvs, xlim = c(0, 1))
# Convergence assessment
ess.beta <- effectiveSize(beta.samples)
summary(ess.beta)
# plot(beta.samples, density = FALSE)
ess.beta.star <- effectiveSize(beta.star.samples)
summary(ess.beta.star)
# plot(kappa.samples, density = FALSE)
ess.kappa <- effectiveSize(kappa.samples)
summary(ess.kappa)
ess.lambda <- effectiveSize(lambda.samples)
summary(ess.lambda[ess.lambda != 0])
# plot(lambda.samples, density = FALSE)

# Devel + apiary ----------------------
load("results/lfMsAbund-devel-apiary-samples.rda")
# GoF assessment
# Plot of true vs. fitted for each species. Looks pretty good for a NB model
for (i in 1:length(sp.names)) {
  plot(c(data.list$y[i, , ]), c(y.rep.means[i, , ]), pch = 19)
  abline(0, 1)
  Sys.sleep(1)
}
# Bayesian p-value for each species
N <- nrow(data.list$y)
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}
hist(bpvs, xlim = c(0, 1))
ess.beta <- effectiveSize(beta.samples)
summary(ess.beta)
# plot(beta.samples, density = FALSE)
ess.beta.star <- effectiveSize(beta.star.samples)
summary(ess.beta.star)
# plot(kappa.samples, density = FALSE)
ess.kappa <- effectiveSize(kappa.samples)
summary(ess.kappa)
ess.lambda <- effectiveSize(lambda.samples)
summary(ess.lambda[ess.lambda != 0])
# plot(lambda.samples, density = FALSE)
