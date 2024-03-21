# 4b-post-process-chains.R: combines the individual chains that were run in parallel
#                           into a single spAbundance model object. This makes it 
#                           easier to assess convergence using the tools that
#                           spAbundance provides.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)
library(coda)
library(abind)

# Load function to combine chians into an spAbundance object --------------
source("code/4a-combineChains.R")

# NOTE: this script will not run as the full model results files are too
#       big for GitHub. If the full results files are desired, please
#       email Jeff Doser (doserjef@msu.edu).

# Directory to save results files -----------------------------------------
# Change as necessary
out.dir <- 'results/'

# Null model --------------------------------------------------------------
load(paste0(out.dir, 'lfMsAbund-null-chain-1-2024-02-04.R'))
out.1 <- out
load(paste0(out.dir, 'lfMsAbund-null-chain-2-2024-02-04.R'))
out.2 <- out
load(paste0(out.dir, 'lfMsAbund-null-chain-3-2024-02-04.R'))
out.3 <- out

out <- combineChains(out.1, out.2, out.3)
save(out, file = paste0(out.dir, 'lfMsAbund-null-all-chains.rda'))

rm(out, out.1, out.2, out.3)
gc()

# Apiary model ------------------------------------------------------------
load(paste0(out.dir, 'lfMsAbund-apiary-chain-1-2024-02-04.R'))
out.1 <- out
load(paste0(out.dir, 'lfMsAbund-apiary-chain-2-2024-02-04.R'))
out.2 <- out
load(paste0(out.dir, 'lfMsAbund-apiary-chain-3-2024-02-04.R'))
out.3 <- out

out <- combineChains(out.1, out.2, out.3)
save(out, file = paste0(out.dir, 'lfMsAbund-apiary-all-chains.rda'))

rm(out, out.1, out.2, out.3)
gc()

# Developed ---------------------------------------------------------------
load(paste0(out.dir, 'lfMsAbund-developed-chain-1-2024-02-04.R'))
out.1 <- out
load(paste0(out.dir, 'lfMsAbund-developed-chain-2-2024-02-04.R'))
out.2 <- out
load(paste0(out.dir, 'lfMsAbund-developed-chain-3-2024-02-04.R'))
out.3 <- out

out <- combineChains(out.1, out.2, out.3)
save(out, file = paste0(out.dir, 'lfMsAbund-developed-all-chains.rda'))

rm(out, out.1, out.2, out.3)
gc()

# Devel + apiary ----------------------------------------------------------
load(paste0(out.dir, 'lfMsAbund-devel-apiary-chain-1-2024-02-04.R'))
out.1 <- out
load(paste0(out.dir, 'lfMsAbund-devel-apiary-chain-2-2024-02-04.R'))
out.2 <- out
load(paste0(out.dir, 'lfMsAbund-devel-apiary-chain-3-2024-02-04.R'))
out.3 <- out

out <- combineChains(out.1, out.2, out.3)
save(out, file = paste0(out.dir, 'lfMsAbund-devel-apiary-all-chains.rda'))

rm(out, out.1, out.2, out.3)
gc()
