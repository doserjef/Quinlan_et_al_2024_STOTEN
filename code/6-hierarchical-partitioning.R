# 6-hierarchical-partitioning.R: script to perform hierarchical partitioning using 
#                                the four candidate model results.
# Authors: Gabriela M. Quinlan and Jeffrey W. Doser
# Note that this code was adapted from code in Zylstra et al. (2021) NEE
# https://github.com/zipkinlab/Zylstra_etal_2021_NEE
library(tidyverse)
library(MCMCvis)

# Load model results ------------------------------------------------------
# Developed + apiary
load("results/lfMsAbund-devel-apiary-samples.rda")
dh <- log.like.by.sp
# Developed only
load('results/lfMsAbund-devel-samples.rda')
d <- log.like.by.sp
# Apiary only
load('results/lfMsAbund-apiary-samples.rda')
h <- log.like.by.sp
# Null model
load('results/lfMsAbund-null-samples.rda')
u <- log.like.by.sp
# Put all log likelihoods by species together in a single matrix
sp.logLik <- rbind(dh, d, h, u)

# Extract species names
sp <- separate(as.data.frame(colnames(beta.samples)[1:length(dh)]), 
               col= "colnames(beta.samples)[1:length(dh)]",
               into = c("blank", "int","sp"))$sp

# Hierarchical partitioning (adapted from Zylstra et al. 2021) ------------
namesOut <- c("dh", "d", "h", "u")

# Final data frame to store all the results
out <- data.frame(matrix(NA, nrow = ncol(sp.logLik), ncol = 2))
colnames(out) <- c("d", "h")

for (j in 1:ncol(sp.logLik)) {
  df <- data.frame(model = namesOut, stringsAsFactors = FALSE)
  df$d <- ifelse(grepl('d',df$model), 1, 0)
  df$h <- ifelse(grepl('h',df$model), 1, 0)
  df$logL <- sp.logLik[,j]
  df[order(-df$logL),]
  
  # Level of hierarchy (0-2 = 0-2 covariate groups, respectively)
  df$hier <- apply(df[,c("d", "h")],1,sum)
  
  # Binary string that indicates which covariates in each model
  df$modelc <- paste0(df$d,df$h)
  
  # Create dataframe with all nested model pairs
  allpairs <- expand.grid(model1 = df$modelc, model0 = df$modelc)
  allpairs <- allpairs[allpairs$model1 != allpairs$model0, ]
  allpairs$hier1 <- df$hier[match(allpairs$model1, df$modelc)]
  allpairs$hier0 <- df$hier[match(allpairs$model0, df$modelc)]
  allpairs <- allpairs[allpairs$hier1 - allpairs$hier0 == 1,]
  allpairs$d.1 <- as.numeric(substr(allpairs$model1, 1, 1))
  allpairs$h.1 <- as.numeric(substr(allpairs$model1, 2, 2))
  
  allpairs$d.0 <- as.numeric(substr(allpairs$model0, 1, 1))
  allpairs$h.0 <- as.numeric(substr(allpairs$model0, 2, 2))
  
  allpairs$d <- allpairs$d.1 + allpairs$d.0
  allpairs$h <- allpairs$h.1 + allpairs$h.0
  
  allpairs$lBoth <- apply(allpairs[, c("d","h")], 1, function(x) sum(x == 2))
  allpairs <- allpairs[allpairs$hier0 == allpairs$lBoth, ]
  
  # Calculate difference in logL for each pair
  allpairs$logL1 <- df$logL[match(allpairs$model1, df$modelc)]
  allpairs$logL0 <- df$logL[match(allpairs$model0, df$modelc)]
  allpairs$diff <- allpairs$logL1 - allpairs$logL0
  
  # Identify covariate group that's different in each pair
  allpairs$param.d <- ifelse(allpairs$d.1 == 1 & allpairs$d.0 == 0, 1, 0)
  allpairs$param.h <- ifelse(allpairs$h.1 == 1 & allpairs$h.0 == 0, 1, 0)
  
  
  # Average logL differences for each covariate group and level of hierarchy
  hpl <- data.frame(expand.grid(param = c("d", "h"),
				hier = unique(allpairs$hier1), stringsAsFactors = FALSE))
  for (i in 1:nrow(hpl)){
    hpl$avg.level[i] <- mean(allpairs$diff[allpairs$hier1 == hpl$hier[i] & 
			     allpairs[, paste0('param.', hpl$param[i])] == 1]) 
  }
  
  # Mean of those averages for each covariate group  
  hp <- data.frame(param=unique(hpl$param))
  for (i in 1:nrow(hp)){
    hp$IC[i] <- mean(hpl$avg.level[hpl$param == hp$param[i]])
    hp$IC.corr[i] <- max(hp$IC[i], 0)
  }
  
  # Relativize values (divide by total)
  hp$IC.perc <- round(hp$IC.corr/sum(hp$IC.corr) * 100,2)
  out[j,"d"]<- hp$IC.perc[1]
  out[j,"h"]<- hp$IC.perc[2]
}

out$sp <- sp 

# Save results ------------------------------------------------------------
write.csv(out, "results/hp-results.csv", row.names = F)

