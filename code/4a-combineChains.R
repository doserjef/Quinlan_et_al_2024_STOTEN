# 4a-combineChains.R: function to process three chains run separately in 
#                     spAbundance into a single object that can then be 
#                     summarized with other spAbundance functions. Note that
#                     this script assumes that you have run three chains. This
#                     is a faster approach to running multiple chains in parallel
#                     as opposed to running them in sequence using the n.chains
#                     argument in spAbundance.
# Author: Jeffrey W. Doser
combineChains <- function(out.1, out.2, out.3) {
  out.full <- list()
  # Combine all MCMC samples together ---------------------------------------
  out.full$n.chains <- 3
  # Community-level effects
  out.full$beta.comm.samples <- mcmc(rbind(out.1$beta.comm.samples, 
  					 out.2$beta.comm.samples, 
  					 out.3$beta.comm.samples))
  # Community-level variances
  out.full$tau.sq.beta.samples <- mcmc(rbind(out.1$tau.sq.beta.samples, 
  					   out.2$tau.sq.beta.samples, 
  					   out.3$tau.sq.beta.samples))
  # Species-level regression coefficients
  out.full$beta.samples <- mcmc(rbind(out.1$beta.samples, 
                                      out.2$beta.samples, 
                                      out.3$beta.samples))
  # Species-specific NB overdispersion parameter (lower = more overdispersion)
  out.full$kappa.samples <- mcmc(rbind(out.1$kappa.samples, 
                                      out.2$kappa.samples, 
                                      out.3$kappa.samples))
  # Overall random effect variances
  out.full$sigma.sq.mu.samples <- mcmc(rbind(out.1$sigma.sq.mu.samples, 
                                      out.2$sigma.sq.mu.samples, 
                                      out.3$sigma.sq.mu.samples))
  # Species-specific random effect values
  out.full$beta.star.samples <- mcmc(rbind(out.1$beta.star.samples, 
                                      out.2$beta.star.samples, 
                                      out.3$beta.star.samples))
  # Factor loadings
  out.full$lambda.samples <- mcmc(rbind(out.1$lambda.samples, 
                                      out.2$lambda.samples, 
                                      out.3$lambda.samples))
  # Factors
  out.full$w.samples <- abind(out.1$w.samples, 
				   out.2$w.samples,
				   out.3$w.samples, along = 1)
  # Fitted values
  out.full$y.rep.samples <- abind(out.1$y.rep.samples, 
                                      out.2$y.rep.samples, 
                                      out.3$y.rep.samples, along = 1)
  # Expected abundance values
  out.full$mu.samples <- abind(out.1$mu.samples, 
                                      out.2$mu.samples, 
                                      out.3$mu.samples, along = 1)
  # Likelihood values (used for WAIC)
  out.full$like.samples <- abind(out.1$like.samples, 
                                      out.2$like.samples, 
                                      out.3$like.samples, along = 1)
  # Get RHat values for the main parameters ---------------------------------
  out.full$rhat <- out.1$rhat
  # beta.comm
  tmp <- mcmc.list(out.1$beta.comm.samples, out.2$beta.comm.samples, out.3$beta.comm.samples)
  out.full$rhat$beta.comm <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
  # tau.sq.beta
  tmp <- mcmc.list(out.1$tau.sq.beta.samples, out.2$tau.sq.beta.samples, out.3$tau.sq.beta.samples)
  out.full$rhat$tau.sq.beta <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
  # beta
  tmp <- mcmc.list(out.1$beta.samples, out.2$beta.samples, out.3$beta.samples)
  out.full$rhat$beta <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
  # kappa
  tmp <- mcmc.list(out.1$kappa.samples, out.2$kappa.samples, out.3$kappa.samples)
  out.full$rhat$kappa <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
  # sigma.sq.mu
  tmp <- mcmc.list(out.1$sigma.sq.mu.samples, out.2$sigma.sq.mu.samples, out.3$sigma.sq.mu.samples)
  out.full$rhat$sigma.sq.mu <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
  # lambda
  tmp <- mcmc.list(out.1$lambda.samples, out.2$lambda.samples, out.3$lambda.samples)
  out.full$rhat$lambda <- as.vector(gelman.diag(tmp, autoburnin = FALSE, 
						multivariate = FALSE)$psrf[, 2])

  
  # Get ESS values for the main parameters ----------------------------------
  out.full$ESS <- list()
  out.full$ESS$beta.comm <- effectiveSize(out.full$beta.comm.samples)
  out.full$ESS$tau.sq.beta <- effectiveSize(out.full$tau.sq.beta.samples)
  out.full$ESS$beta <- effectiveSize(out.full$beta.samples)
  out.full$ESS$kappa <- effectiveSize(out.full$kappa.samples)
  out.full$ESS$sigma.sq.mu <- effectiveSize(out.full$sigma.sq.mu.samples)
  out.full$ESS$lambda <- effectiveSize(out.full$lambda.samples)
  
  # Other stuff -------------------------------------------------------------
  # This stuff is used under the hood for use with summaries, plotting, 
  # prediction, WAIC,  PPCs. 
  # Names of random effects, which are needed to keep track of stuff if 
  # doing any prediction
  out.full$re.level.names <- out.1$re.level.names
  # Design matrix for fixed effects
  out.full$X <- out.1$X
  # Design matrix for random effects
  out.full$X.re <- out.1$X.re
  # The count data
  out.full$y <- out.1$y
  # Information on the call to the function
  out.full$call <- out.1$call
  # Overall number of samples run per chain
  out.full$n.samples <- out.1$n.samples
  # Names of columns in C
  out.full$x.names <- out.1$x.names
  # Species names
  out.full$sp.names <- out.1$sp.names
  # Number of posterior samples saved for each chain
  out.full$n.post <- out.1$n.post
  # Thinning rate
  out.full$n.thin <- out.1$n.thin
  # Amount of burn-in
  out.full$n.burn <- out.1$n.burn
  # Number of chains
  out.full$n.chains <- 3
  # Distribution used for abundance (NB or Poisson)
  out.full$dist <- out.1$dist
  # Names of columns with random effects
  out.full$re.cols <- out.1$re.cols
  # Logical indicating if there were any random effects
  out.full$muRE <- out.1$muRE
  # Offset (needed for posterior predictive checks)
  out.full$offset <- matrix(1, ncol(out.full$y), dim(out.full$y)[3])
  # Number of factors
  out.full$q <- out.1$q
  # Run time
  # Setting the "overall" run time to be the longest run time across the three chains
  tmp <- which.max(c(out.1$run.time[3], out.2$run.time[3], out.3$run.time[3]))
  if (tmp == 1) out.full$run.time <- out.1$run.time
  if (tmp == 2) out.full$run.time <- out.2$run.time
  if (tmp == 3) out.full$run.time <- out.3$run.time
 
  # Assign class 
  class(out.full) <- class(out.1) 

  return(out.full)
} 
