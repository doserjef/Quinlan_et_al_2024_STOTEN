# 3c-main-apiary.R: script to run a JSDM using spAbundance using apiary 
#                   density as a covariate in the model. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# Get chain number from command line run ----------------------------------
# Uncomment this line if running scripts from command line
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# Alternatively, if not running the script from the command line, you can 
# manually change the script number below to run different chains:
# chain <- 1
# Or, can use the n.chains function in spOccupancy (for sequential runs of
# chains).
if(length(chain) == 0) base::stop('Need to tell spAbundance the chain number')

# Read in data ------------------------------------------------------------
load('data/spAbundance-data.rda')

# Prep model for spOccupancy ----------------------------------------------
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   sigma.sq.mu.ig = list(a = 0.1, b = 0.1))

# MCMC stuff
n.batch <- 5000
batch.length <- 25
n.burn <- 65000
n.thin <- 20
n.chains <- 1

# Load initial values
load("results/inits-apiary.rda")
inits.list$w <- NULL
inits.list$lambda <- NULL

# Run the model -----------------------------------------------------------
out <- lfMsAbund(formula = ~ scale(apiary.density) + (1 | day.random) + 
                             factor(year) + factor(protocol) + effort:protocol,
                 data = data.list, 
                 n.batch = n.batch, 
                 batch.length = batch.length, 
                 inits = inits.list,
                 priors = prior.list, 
                 n.factors = 8,
                 family = 'NB',
                 verbose = TRUE, 
                 n.report = 50,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = n.chains) 

summary(out)

# Save to hard drive ------------------------------------------------------
# NOTE: this file is too large for GitHub.
save(out, file = paste0("/mnt/disk4/jeff/QDKG23/results/lfMsAbund-apiary-chain-", 
                        chain, "-", Sys.Date(), ".R"))
