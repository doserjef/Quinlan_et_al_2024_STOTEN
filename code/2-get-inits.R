# 2-get-inits.R: script to extract initial values from a previous model fit 
#                to aid in convergence.
#                NOTE: previous model fit results files are not available on GitHub
#                      so this script will not run, but it shows how the initial 
#                      values are generated, which we use to aid in convergence.
# Author: Jeffrey W. Doser
rm(list = ls())
library(coda)
library(spAbundance)

# Apiary only -------------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-apiary-2023-12-10.R")

n.sp <- nrow(data.list$y)
p.abund <- dim(out$X)[3]
q <- out$q
J <- dim(out$X)[1]
beta.comm.inits <- apply(out$beta.comm.samples[1:3000, ], 2, median)
tau.sq.beta.inits <- apply(out$tau.sq.beta.samples[1:3000, ], 2, median)
beta.inits <- matrix(apply(out$beta.samples[1:3000, ], 2, median), n.sp, p.abund)
sigma.sq.mu.inits <- median(out$sigma.sq.mu.samples[1:3000, ])
kappa.inits <- apply(out$kappa.samples[1:3000, ], 2, median)
lambda.inits <- matrix(apply(out$lambda.samples[1:3000, ], 2, median), n.sp, q)
w.inits <- apply(out$w.samples[1:3000, , ], c(2, 3), median)
inits.list <- list(beta.comm = beta.comm.inits, tau.sq.beta = tau.sq.beta.inits,
		   beta = beta.inits, sigma.sq.mu = sigma.sq.mu.inits, 
		   kappa = kappa.inits, lambda = lambda.inits, w = w.inits)
save(inits.list, file = '/mnt/disk4/jeff/QDKG23/results/inits-apiary.rda')

# Developed only ----------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-developed-2023-12-10.R")
n.sp <- nrow(data.list$y)
p.abund <- dim(out$X)[3]
q <- out$q
J <- dim(out$X)[1]
beta.comm.inits <- apply(out$beta.comm.samples[1:3000, ], 2, median)
tau.sq.beta.inits <- apply(out$tau.sq.beta.samples[1:3000, ], 2, median)
beta.inits <- matrix(apply(out$beta.samples[1:3000, ], 2, median), n.sp, p.abund)
sigma.sq.mu.inits <- median(out$sigma.sq.mu.samples[1:3000, ])
kappa.inits <- apply(out$kappa.samples[1:3000, ], 2, median)
lambda.inits <- matrix(apply(out$lambda.samples[1:3000, ], 2, median), n.sp, q)
w.inits <- apply(out$w.samples[1:3000, , ], c(2, 3), median)
inits.list <- list(beta.comm = beta.comm.inits, tau.sq.beta = tau.sq.beta.inits,
		   beta = beta.inits, sigma.sq.mu = sigma.sq.mu.inits, 
		   kappa = kappa.inits, lambda = lambda.inits, w = w.inits)
save(inits.list, file = '/mnt/disk4/jeff/QDKG23/results/inits-developed.rda')

# Developed + apiary ------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-devel-apiary-2023-12-10.R")

n.sp <- nrow(data.list$y)
p.abund <- dim(out$X)[3]
q <- out$q
J <- dim(out$X)[1]
beta.comm.inits <- apply(out$beta.comm.samples[1:3000, ], 2, median)
tau.sq.beta.inits <- apply(out$tau.sq.beta.samples[1:3000, ], 2, median)
beta.inits <- matrix(apply(out$beta.samples[1:3000, ], 2, median), n.sp, p.abund)
sigma.sq.mu.inits <- median(out$sigma.sq.mu.samples[1:3000, ])
kappa.inits <- apply(out$kappa.samples[1:3000, ], 2, median)
lambda.inits <- matrix(apply(out$lambda.samples[1:3000, ], 2, median), n.sp, q)
w.inits <- apply(out$w.samples[1:3000, , ], c(2, 3), median)
inits.list <- list(beta.comm = beta.comm.inits, tau.sq.beta = tau.sq.beta.inits,
		   beta = beta.inits, sigma.sq.mu = sigma.sq.mu.inits, 
		   kappa = kappa.inits, lambda = lambda.inits, w = w.inits)
save(inits.list, file = '/mnt/disk4/jeff/QDKG23/results/inits-devel-apiary.rda')

# Null model --------------------------------------------------------------
load("data/spAbundance-data.rda")
load("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-null-2023-12-10.R")

n.sp <- nrow(data.list$y)
p.abund <- dim(out$X)[3]
q <- out$q
J <- dim(out$X)[1]
beta.comm.inits <- apply(out$beta.comm.samples[1:3000, ], 2, median)
tau.sq.beta.inits <- apply(out$tau.sq.beta.samples[1:3000, ], 2, median)
beta.inits <- matrix(apply(out$beta.samples[1:3000, ], 2, median), n.sp, p.abund)
sigma.sq.mu.inits <- median(out$sigma.sq.mu.samples[1:3000, ])
kappa.inits <- apply(out$kappa.samples[1:3000, ], 2, median)
lambda.inits <- matrix(apply(out$lambda.samples[1:3000, ], 2, median), n.sp, q)
w.inits <- apply(out$w.samples[1:3000, , ], c(2, 3), median)
inits.list <- list(beta.comm = beta.comm.inits, tau.sq.beta = tau.sq.beta.inits,
		   beta = beta.inits, sigma.sq.mu = sigma.sq.mu.inits, 
		   kappa = kappa.inits, lambda = lambda.inits, w = w.inits)
save(inits.list, file = '/mnt/disk4/jeff/QDKG23/results/inits-null.rda')
