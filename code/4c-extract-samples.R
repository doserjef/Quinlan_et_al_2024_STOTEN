# 4c-extract-samples.R: script to extract the relevant information from the 
#                       full spAbundance model objects and save in a smaller, 
#                       more manageable file. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(coda)
library(spAbundance)

# NOTE: this script will not run as the full results files are too large for
#       GitHub. If full results files are desired, please email Jeff Doser
#       (doserjef@msu.edu). If you have run all preceding scripts, then this 
#       script should run successfully. 
# Null model --------------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-null-all-chains.rda")

rhat.vals <- unlist(out$rhat)
beta.samples <- out$beta.samples
beta.star.samples <- out$beta.star.samples
kappa.samples <- out$kappa.samples
# Mean fitted values
y.rep.means <- apply(out$y.rep.samples, c(2, 3, 4), mean)
lambda.samples <- out$lambda.samples
waic.model <- waicAbund(out, by.sp = TRUE)

# Run posterior predictive check
ppc.out <- ppcAbund(out, fit.stat = 'freeman-tukey', group = 0)

fit.y <- ppc.out$fit.y
fit.y.rep <- ppc.out$fit.y.rep
N <- ncol(fit.y)
log.like.by.sp <- rep(NA, N)
for (i in 1:N) {
  print(paste0('Species ', i, ' out of ', N))
  log.like.by.sp[i] <- mean(apply(log(out$like.samples[, i, , ]), 1, sum, na.rm = TRUE))
}
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}

save(beta.samples, beta.star.samples, kappa.samples, rhat.vals,
     y.rep.means, lambda.samples, bpvs, waic.model, log.like.by.sp,
     file = '/mnt/disk4/jeff/QDKG23/results/lfMsAbund-null-samples.rda')

rm(y.rep.means, ppc.out, fit.y, fit.y.rep, out)
gc()

# Apiary only -------------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-apiary-all-chains.rda")

rhat.vals <- unlist(out$rhat)
sum(rhat.vals > 1.1)
rhat.vals[which(rhat.vals > 1.1)]
beta.samples <- out$beta.samples
beta.star.samples <- out$beta.star.samples
kappa.samples <- out$kappa.samples
# Mean fitted values
y.rep.means <- apply(out$y.rep.samples, c(2, 3, 4), mean)
lambda.samples <- out$lambda.samples
waic.model <- waicAbund(out, by.sp = TRUE)

# Run posterior predictive check
ppc.out <- ppcAbund(out, fit.stat = 'freeman-tukey', group = 0)

fit.y <- ppc.out$fit.y
fit.y.rep <- ppc.out$fit.y.rep
N <- nrow(data.list$y)
log.like.by.sp <- rep(NA, N)
for (i in 1:N) {
  print(paste0('Species ', i, ' out of ', N))
  log.like.by.sp[i] <- mean(apply(log(out$like.samples[, i, , ]), 1, sum, na.rm = TRUE))
}
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}

save(beta.samples, beta.star.samples, kappa.samples, rhat.vals,
     y.rep.means, lambda.samples, bpvs, waic.model, log.like.by.sp,
     file = '/mnt/disk4/jeff/QDKG23/results/lfMsAbund-apiary-samples.rda')
rm(y.rep.means, ppc.out, fit.y, fit.y.rep, out)
gc()
# Developed only ----------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-developed-all-chains.rda")

rhat.vals <- unlist(out$rhat)
sum(rhat.vals > 1.1)
rhat.vals[which(rhat.vals > 1.1)]
beta.samples <- out$beta.samples
beta.star.samples <- out$beta.star.samples
kappa.samples <- out$kappa.samples
# Mean fitted values
y.rep.means <- apply(out$y.rep.samples, c(2, 3, 4), mean)
lambda.samples <- out$lambda.samples
waic.model <- waicAbund(out, by.sp = TRUE)

# Run posterior predictive check
ppc.out <- ppcAbund(out, fit.stat = 'freeman-tukey', group = 0)

fit.y <- ppc.out$fit.y
fit.y.rep <- ppc.out$fit.y.rep
N <- nrow(data.list$y)
log.like.by.sp <- rep(NA, N)
for (i in 1:N) {
  print(paste0('Species ', i, ' out of ', N))
  log.like.by.sp[i] <- mean(apply(log(out$like.samples[, i, , ]), 1, sum, na.rm = TRUE))
}
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}

save(beta.samples, beta.star.samples, kappa.samples, rhat.vals,
     y.rep.means, lambda.samples, bpvs, waic.model, log.like.by.sp,
     file = '/mnt/disk4/jeff/QDKG23/results/lfMsAbund-devel-samples.rda')
rm(y.rep.means, ppc.out, fit.y, fit.y.rep, out)
gc()

# Developed + apiary ------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-devel-apiary-all-chains.rda")

rhat.vals <- unlist(out$rhat)
sum(rhat.vals > 1.1)
rhat.vals[which(rhat.vals > 1.1)]
beta.samples <- out$beta.samples
beta.star.samples <- out$beta.star.samples
kappa.samples <- out$kappa.samples
# Mean fitted values
y.rep.means <- apply(out$y.rep.samples, c(2, 3, 4), mean)
lambda.samples <- out$lambda.samples
waic.model <- waicAbund(out, by.sp = TRUE)

# Run posterior predictive check
ppc.out <- ppcAbund(out, fit.stat = 'freeman-tukey', group = 0)

fit.y <- ppc.out$fit.y
fit.y.rep <- ppc.out$fit.y.rep
N <- nrow(data.list$y)
log.like.by.sp <- rep(NA, N)
for (i in 1:N) {
  print(paste0('Species ', i, ' out of ', N))
  log.like.by.sp[i] <- mean(apply(log(out$like.samples[, i, , ]), 1, sum, na.rm = TRUE))
}
bpvs <- rep(NA, N)
for (i in 1:N) {
 bpvs[i] <- mean(fit.y.rep[, i] > fit.y[, i])
}

save(beta.samples, beta.star.samples, kappa.samples, rhat.vals,
     y.rep.means, lambda.samples, bpvs, waic.model, log.like.by.sp, 
     file = '/mnt/disk4/jeff/QDKG23/results/lfMsAbund-devel-apiary-samples.rda')
rm(y.rep.means, ppc.out, fit.y, fit.y.rep, out)
gc()

