# 7a-summary.R: script to summarize all analysis results and generate all figures
#               provided in the manuscript.
# Authors: Jeffrey W. Doser and Gabriela M. Quinlan
rm(list = ls())
library(tidyverse)
library(spAbundance)
library(coda)
library(MCMCvis)
library(ggpubr)
library(sf)
library(fields)
library(patchwork)
library(stars)

# NOTE: the results files need to be generated prior to running this script.
# Load data ---------------------------------------------------------------
load("data/spAbundance-data.rda")
sp.names <- dimnames(data.list$y)[[1]]

# Load model results from devel only and apiary only models ---------------
# Apiary only model
load("results/lfMsAbund-apiary-samples.rda")
beta.apiary.samples <- MCMCchains(beta.samples, params = 'apiary.density', exact = FALSE)
beta.apiary.prob.pos <- apply(beta.apiary.samples, 2, function(a) mean(a > 0))
N <- length(sp.names)
beta.pan.net.samples <- MCMCchains(beta.samples, params = 'panTrap', exact = FALSE)[, 1:N]
beta.vane.net.samples <- MCMCchains(beta.samples, params = 'vaneTrap', exact = FALSE)[, 1:N]
# Notice the negative in the function is to get this as the prob that hand trap
# (which is the reference level) is higher than pan trap (or vane trap)
beta.hand.net.greater.than.pan.prob.pos <- apply(beta.pan.net.samples, 2, function(a) mean(a < 0))
beta.hand.net.greater.than.vane.prob.pos <- apply(beta.vane.net.samples, 2, function(a) mean(a < 0))
# Developed only model
load("results/lfMsAbund-devel-samples.rda")
beta.devel.samples <- MCMCchains(beta.samples, params = 'dev', exact = FALSE)
beta.devel.prob.pos <- apply(beta.devel.samples, 2, function(a) mean(a > 0))

# Save certain results described in the paper in a csv
probs.df <- data.frame(genus = sp.names, 
                       apiary.prob.pos = beta.apiary.prob.pos,
                       hand.greater.than.pan.prob.pos = beta.hand.net.greater.than.pan.prob.pos,
                       hand.greater.than.vane.prob.pos = beta.hand.net.greater.than.vane.prob.pos,
                       devel.prob.pos = beta.devel.prob.pos)
rownames(probs.df) <- NULL
write.csv(probs.df, file = 'results/beta-subset-prob-results.csv', row.names = FALSE)

# Generate caterpillar plots ----------------------------------------------
# Apiary effect -----------------------
beta.apiary.df <- data.frame(med = apply(beta.apiary.samples, 2, median), 
                             low = apply(beta.apiary.samples, 2, quantile, 0.25),
                             high = apply(beta.apiary.samples, 2, quantile, 0.75), 
                             lowest = apply(beta.apiary.samples, 2, quantile, 0.025),
                             highest = apply(beta.apiary.samples, 2, quantile, 0.975), 
                             prob.neg = 1 - beta.apiary.prob.pos,
                             raw.counts = apply(data.list$y, 1, sum, na.rm = TRUE),
                             sp = sp.names)

plot.order <- sp.names[order(beta.apiary.df$med)]
apiary.cat.plot <- beta.apiary.df %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = med, y = sp, fill = prob.neg)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = lowest, y = sp, xend = highest, yend = sp),
                 lineend = 'round', linewidth = 1, col = 'lightgray') +
    geom_segment(aes(x = low, y = sp, xend = high, yend = sp),
                 lineend = 'round', linewidth = 2, col = 'darkgray') +
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, high = '#B2182B', mid = 'white', low = '#2166AC',
                         na.value = NA, limits = c(0, 1)) +
    theme_bw(base_size = 17) +
    labs(x = 'Apiary Density Effect Size',
         y = 'Genus', fill = 'P(effect < 0)', title = '(A)') +
    theme(text = element_text(family = 'LM Roman 10'), 
          panel.grid = element_blank(), 
          legend.position.inside = c(0.85, 0.25), 
          legend.position = 'inside',
          plot.title = element_text(size = 17),
          legend.background = element_rect(fill = NA))
apiary.cat.plot
# Devel effect ------------------------
beta.devel.df <- data.frame(med = apply(beta.devel.samples, 2, median), 
                            low = apply(beta.devel.samples, 2, quantile, 0.25),
                            high = apply(beta.devel.samples, 2, quantile, 0.75), 
                            lowest = apply(beta.devel.samples, 2, quantile, 0.025),
                            highest = apply(beta.devel.samples, 2, quantile, 0.975), 
                            prob.neg = 1 - beta.devel.prob.pos,
                            sp = sp.names)
plot.order <- sp.names[order(beta.devel.df$med)]
devel.cat.plot <- beta.devel.df %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = med, y = sp, fill = prob.neg)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = lowest, y = sp, xend = highest, yend = sp),
                 lineend = 'round', linewidth = 1, col = 'lightgray') +
    geom_segment(aes(x = low, y = sp, xend = high, yend = sp),
                 lineend = 'round', linewidth = 2, col = 'darkgray') +
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, high = '#B2182B', mid = 'white', low = '#2166AC',
                         na.value = NA, limits = c(0, 1)) +
    theme_bw(base_size = 17) +
    labs(x = 'Developed Land Effect Size',
         y = 'Genus', fill = 'P(effect < 0)', title = '(B)') +
    theme(text = element_text(family = 'LM Roman 10'), 
          panel.grid = element_blank(), 
          legend.position.inside = c(0.85, 0.25), 
          legend.position = 'inside',
          plot.title = element_text(size = 17),
          legend.background = element_rect(fill = NA))
devel.cat.plot

fig.2.plot <- ggarrange(apiary.cat.plot, devel.cat.plot, nrow = 1)
ggsave(fig.2.plot, file = 'figures/Figure-2.png', device = 'png', 
       units = 'in', width = 14, height = 9)

# More in-depth effects of apiary density ---------------------------------
load('results/lfMsAbund-apiary-samples.rda')
# Apiary density posterior samples
beta.apiary.samples <- MCMCchains(beta.samples, params = 'apiary.density', exact = FALSE)
# Intecept posterior samples
beta.int.samples <- MCMCchains(beta.samples, params = 'Intercept', exact = FALSE)
# Pan trap effect posterior samples
beta.pan.samples <- beta.pan.net.samples
# Pan trap effect means
beta.pan.means <- apply(beta.pan.samples, 2, mean)
# Vane trap effect posterior samples
beta.vane.samples <- beta.vane.net.samples
# Vane trap effect means
beta.vane.means <- apply(beta.vane.samples, 2, mean)
# A vector that indicates which sampling method was the best for each individual genera
method.ind <- vector(mode = 'numeric', length = N)
for (i in 1:N) {
  if (beta.pan.means[i] > beta.vane.means[i] & beta.pan.means[i] > 0) {
    method.ind[i] <- 'pan'
  }
  if (beta.vane.means[i] > beta.pan.means[i] & beta.vane.means[i] > 0) {
    method.ind[i] <- 'vane'
  }
  if (beta.vane.means[i] < 0 & beta.pan.means[i] < 0) {
    method.ind[i] <- 'hand'
  }
}
# A large vector of potential apiary values to predict over
n.pred <- 1000
apiary.pred.vals <- seq(min(data.list$covs$apiary.density), 
		   max(data.list$covs$apiary.density), 
		   length.out = n.pred)
# Scale by values used to fit model
apiary.pred.vals.s <- (apiary.pred.vals - mean(data.list$covs$apiary.density)) / 
                      sd(data.list$covs$apiary.density)
# Only predict across the smallest and largest apiary density value
pred.vals <- c(apiary.pred.vals[1], apiary.pred.vals[n.pred])
pred.vals.s <- c(apiary.pred.vals.s[1], apiary.pred.vals.s[n.pred])
n.post <- nrow(beta.apiary.samples)

# Generate posterior samples for predicted abundances for each genus
# at the highest and lowest values
pred.vals.samples <- array(NA, dim = c(n.post, length(pred.vals), N))
for (i in 1:N) {
  for (j in 1:length(pred.vals)) {
    if (method.ind[i] == 'hand') {
      pred.vals.samples[, j, i] <- exp(beta.int.samples[, i] + beta.apiary.samples[, i] * pred.vals.s[j])
    }
    if (method.ind[i] == 'pan') {
      pred.vals.samples[, j, i] <- exp(beta.int.samples[, i] + 
				       beta.apiary.samples[, i] * pred.vals.s[j] +
                                       beta.pan.samples[, i])
    }
    if (method.ind[i] == 'vane') {
      pred.vals.samples[, j, i] <- exp(beta.int.samples[, i] + 
				       beta.apiary.samples[, i] * pred.vals.s[j] +
                                       beta.vane.samples[, i])
    }
  }
}
# Median values
pred.vals.meds <- apply(pred.vals.samples, c(2, 3), median)
# Lower bound of 95% CI
pred.vals.lowest <- apply(pred.vals.samples, c(2, 3), quantile, 0.025)
# Lower bound of 75% CI
pred.vals.low <- apply(pred.vals.samples, c(2, 3), quantile, 0.25)
# Upper bound of 95% CI
pred.vals.highest <- apply(pred.vals.samples, c(2, 3), quantile, 0.975)
# Upper bound of 75% CI
pred.vals.high <- apply(pred.vals.samples, c(2, 3), quantile, 0.75)

# Percent change in medians from highest to lowest
percent.change <- pred.vals.meds[2, ] / pred.vals.meds[1, ] * 100
pc.df <- data.frame(change = percent.change, 
                    sp = sp.names)
pc.df %>%
  filter(sp %in% c('Svastra', 'Melitoma', 'Triepeolus', 'Augochloropsis')) %>%
  mutate(change = 100 - change)


plot.df <- data.frame(med = c(pred.vals.meds),
                      low = c(pred.vals.low),
                      lowest = c(pred.vals.lowest),
                      high = c(pred.vals.high),
                      highest = c(pred.vals.highest),
                      type = factor(rep(c('Lowest', 'Highest'), times = N), 
                                    levels = c('Lowest', 'Highest')),
                      sp = rep(sp.names, each = nrow(pred.vals.meds)))
# Save plot to csv for summarizing in table
write.csv(plot.df, file = 'results/effect-sizes-apiary-pred.csv', row.names = FALSE)
# Subset to 4 genera with decent support for effect of apiary density
small.plot.df <- plot.df %>% 
  filter(sp %in% c('Svastra', 'Melitoma', 'Triepeolus', 'Augochloropsis'))
small.plot.df$sp[which(small.plot.df$sp == 'Svastra')] <- '(A) Svastra'
small.plot.df$sp[which(small.plot.df$sp == 'Melitoma')] <- '(B) Melitoma'
small.plot.df$sp[which(small.plot.df$sp == 'Triepeolus')] <- '(C) Triepeolus'
small.plot.df$sp[which(small.plot.df$sp == 'Augochloropsis')] <- '(D) Augochloropsis'
# Generate plot for 4 genera with decent support for effect of apiary density
fig.3.plot <- small.plot.df %>%
    ggplot(aes(x = type, y = med)) +
      geom_segment(aes(x = type, y = lowest, xend = type, yend = highest),
          	 lineend = 'round', linewidth = 1, col = 'lightgray') +
      geom_segment(aes(x = type, y = low, xend = type, yend = high),
                   lineend = 'round', linewidth = 2, col = 'darkgrey') +
      geom_point() +
      geom_line(aes(group = sp)) +
      facet_wrap(vars(sp), scales = 'free_y') +
      theme_bw(base_size = 16) +
      labs(x = 'Apiary density', y = 'Expected relative abundance') + 
      theme(text = element_text(family = 'LM Roman 10'), 
            panel.grid = element_blank(), 
            legend.position.inside = c(0.85, 0.25), 
            legend.position = 'inside',
            plot.title = element_text(size = 17),
            legend.background = element_rect(fill = NA))
ggsave(plot = fig.3.plot, file = 'figures/Figure-3.png', device = 'png', 
       units = 'in', height = 6, width = 8)

# Generate plot for all genera
fig.s4.plot <- plot.df %>%
    ggplot(aes(x = type, y = med)) +
      geom_segment(aes(x = type, y = lowest, xend = type, yend = highest),
          	 lineend = 'round', linewidth = 1, col = 'lightgray') +
      geom_segment(aes(x = type, y = low, xend = type, yend = high),
                   lineend = 'round', linewidth = 2, col = 'darkgray') +
      geom_point() +
      geom_line(aes(group = sp)) +
      facet_wrap(vars(sp), scales = 'free_y') +
      theme_bw() +
      labs(x = 'Apiary density', y = 'Expected relative abundance') + 
      theme(text = element_text(family = 'LM Roman 10'), 
            panel.grid = element_blank(), 
            legend.position = c(0.85, 0.25), 
            plot.title = element_text(size = 17),
            legend.background = element_rect(fill = NA))
ggsave(plot = fig.s4.plot, file = 'figures/Figure-S4.png', device = 'png', 
       units = 'in', height = 8, width = 12)

# Visualize hierarchical partitioning results -----------------------------
# Load hierarchical partitioning results
hp.out <- read.csv("results/hp-results.csv")
sel.order <- hp.out %>%
  arrange(h)

plot.df <- hp.out %>%
  rename(species = sp) %>%
  pivot_longer(cols = d:h, values_to = "ic", names_to = "var") %>%
  mutate(labels = factor(species, levels = sel.order$sp, ordered = TRUE)) %>%
  na.omit() %>%
  mutate(fill = ifelse(var == "h", "No Support (HB Density)", "No Support (Developed)")) %>%
  # negHB species are species with .8 prob of negative apiary effect
  mutate(fill = ifelse(fill == "No Support (HB Density)" & 
                               species %in% sp.names[which(beta.apiary.prob.pos < .2)], 
                       "Negative (HB Density)", fill)) %>%
  # posHB species are speices with .8 prob of positive apiary effect
  mutate(fill = ifelse(fill == "No Support (HB Density)" & 
                               species %in% sp.names[which(beta.apiary.prob.pos > .8)], 
                       "Positive (HB Density)", fill)) %>%   
  # negDev species are species with .8 prob of negative developed effect
  mutate(fill = ifelse(fill == "No Support (Developed)" & 
                               species %in% sp.names[which(beta.devel.prob.pos < .2)], 
                       "Negative (Developed)", fill)) %>%    
  # posDev species are species with .8 prob of positive developed effect
  mutate(fill = ifelse(fill == "No Support (Developed)" & 
                               species %in% sp.names[which(beta.devel.prob.pos > .8)], 
                       "Positive (Developed)", fill)) %>%
  mutate(fill = factor(fill, levels = c("No Support (Developed)", "Negative (Developed)", 
                                        "Positive (Developed)", 
                                        "No Support (HB Density)", "Negative (HB Density)", 
                                        "Positive (HB Density)")))

hp.plot <- ggplot(data = plot.df, aes(x = labels, y = ic)) +
  geom_bar(stat = "identity", aes(fill = fill)) +
  coord_flip() +
  geom_hline(yintercept= 50, lty = 2) +
  scale_fill_manual(values = c("grey30", "darkgreen", "lightgreen", "grey", "orange", "gold"))+
  theme_bw(base_size = 14) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(y = 'Relative Importance (%)', x = 'Genus', fill = '') + 
  theme(text = element_text(family = 'LM Roman 10'), 
        axis.ticks.y = element_blank())
hp.plot
ggsave(hp.plot, file = 'figures/Figure-4.png', device = 'png', 
       units = 'in', width = 7, height = 7)

# Visualize sampling protocol differences across genera -------------------
# Pan trapping relative to handnet
beta.pan.df <- data.frame(med = apply(beta.pan.net.samples, 2, median), 
                          low = apply(beta.pan.net.samples, 2, quantile, 0.25),
                          high = apply(beta.pan.net.samples, 2, quantile, 0.75), 
                          lowest = apply(beta.pan.net.samples, 2, quantile, 0.025),
                          highest = apply(beta.pan.net.samples, 2, quantile, 0.975), 
                          prob.neg = beta.hand.net.greater.than.pan.prob.pos,
                          sp = sp.names)
plot.order <- sp.names[order(beta.pan.df$med)]
pan.plot <- beta.pan.df %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = med, y = sp, fill = prob.neg)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = lowest, y = sp, xend = highest, yend = sp),
                 lineend = 'round', linewidth = 1, col = 'lightgray') +
    geom_segment(aes(x = low, y = sp, xend = high, yend = sp),
                 lineend = 'round', linewidth = 2, col = 'darkgray') +
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, high = '#B2182B', mid = 'white', low = '#2166AC',
                         na.value = NA, limits = c(0, 1)) +
    theme_bw(base_size = 17) +
    labs(x = 'Effect of pan trapping relative to hand netting',
         y = 'Genus', fill = 'P(effect < 0)', title = '(A)') +
    theme(text = element_text(family = 'LM Roman 10'), 
          panel.grid = element_blank(), 
          legend.position.inside = c(0.85, 0.25), 
          legend.position = 'inside',
          plot.title = element_text(size = 17),
          legend.background = element_rect(fill = NA))
pan.plot

# Vane trapping relative to hand netting
beta.vane.df <- data.frame(med = apply(beta.vane.net.samples, 2, median), 
                           low = apply(beta.vane.net.samples, 2, quantile, 0.25),
                           high = apply(beta.vane.net.samples, 2, quantile, 0.75), 
                           lowest = apply(beta.vane.net.samples, 2, quantile, 0.025),
                           highest = apply(beta.vane.net.samples, 2, quantile, 0.975), 
                           prob.neg = beta.hand.net.greater.than.vane.prob.pos,
                           sp = sp.names)
plot.order <- sp.names[order(beta.vane.df$med)]
vane.plot <- beta.vane.df %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = med, y = sp, fill = prob.neg)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = lowest, y = sp, xend = highest, yend = sp),
                 lineend = 'round', linewidth = 1, col = 'lightgray') +
    geom_segment(aes(x = low, y = sp, xend = high, yend = sp),
                 lineend = 'round', linewidth = 2, col = 'darkgray') +
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, high = '#B2182B', mid = 'white', low = '#2166AC',
                         na.value = NA, limits = c(0, 1)) +
    theme_bw(base_size = 17) +
    labs(x = 'Effect of vane trapping relative to hand netting',
         y = 'Genus', fill = 'P(effect < 0)', title = '(B)') +
    theme(text = element_text(family = 'LM Roman 10'), 
          panel.grid = element_blank(), 
          legend.position = c(0.85, 0.25), 
          plot.title = element_text(size = 17),
          legend.background = element_rect(fill = NA))
vane.plot

fig.4.plot <- ggarrange(pan.plot, vane.plot, nrow = 1)
ggsave(fig.4.plot, file = 'figures/Figure-5.png', device = 'png', 
       units = 'in', width = 14, height = 9)

# Plot of apiary density metric (Figure S1) -------------------------------
# Plot of the function ----------------
n <- 10000
x <- seq(0, 6, length.out = n)
y <- 1 / exp(x / 2)
plot.df <- data.frame(y = y, x = x)
inv.exp.plot <- ggplot(plot.df, aes(x = x, y = y)) +
  geom_line(linewidth = 0.8, lineend = 'round') +
  theme_bw(base_size = 16) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 6)) +
  labs(x = 'Distance from Apiary(km)', y = 'Value', title = '(A)') +
  theme(text = element_text(family = 'LM Roman 10'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16))
# Plot of a single point in space -----
my.point <- st_point(c(0, 0))
point.buffered <- st_buffer(my.point, 6)
# Grid up the circle
grid.pred <- st_as_stars(st_bbox(point.buffered), dx = .1, dy = .1)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred,
                           coords = c('x', 'y'))

# Intersect with region of interest
# Takes a few minutes
coords.pred.sf <- st_intersection(coords.pred.sf, point.buffered)
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))
dist <- rep(NA, nrow(coords.0))
for (j in 1:nrow(coords.0)) {
  dist[j] <- st_distance(my.point, coords.pred.sf[j, ])
}
plot.apiary.density <- 1 / exp(dist / 2)

df <- data.frame(density = plot.apiary.density,
                 dist = dist,
                 X = coords.0[, 1],
                 Y = coords.0[, 2])
example.apiary.plot <- ggplot(df, aes(X, Y, fill = density)) +
  geom_raster() + 
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(fill = "Value", x = 'Easting', y = 'Northing', title = '(B)') +
  annotate("path", x = 2 * cos(seq(0, 2 * pi, length.out = 100)),
           y = 2 * sin(seq(0, 2 * pi, length.out = 100)), col = 'white') +
  theme_bw(base_size = 16) +
  theme(legend.position = "right",
        text = element_text(family = 'LM Roman 10'),
        plot.title = element_text(size = 16),
        panel.grid = element_blank())

# Generate plot of two different apiary density metrics -------------------
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs" 
coords.aea <- st_as_sf(as.data.frame(data.list$coords), 
                       coords = c('X', 'Y'), 
                       crs = my.crs)
# Compare current apiary density metric to a static 3.5km metric.
# Get ccoordinates in AEA
coords.buffered <- st_buffer(coords.aea, dist = 3.5)
maryland <- st_as_sf(maps::map('state', region = c('maryland'), 
                               fill = TRUE, plot = FALSE))
maryland.wgs <- maryland %>%
  st_transform(4326)
# Check to make sure the complete circle is within Maryland boundaries
maryland.aea <- st_transform(maryland.wgs, my.crs)
# Load apiary data --------------------------------------------------------
# NOTE: this code won't run given the propietary nature of the apiary locations.
apiary <- read.csv("data/geocoded_apiary_locations1718.csv")
apiary <- apiary %>%
  filter(lat != 'ERROR', lon != 'ERROR', ST == 'MD')
coords.apiary <- apiary %>%
  dplyr::select(lon, lat) %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon))
coords.apiary.sf <- st_as_sf(coords.apiary, 
                             coords = c('lon', 'lat'), 
                             crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
apiary.aea <- coords.apiary.sf %>%
  st_transform(my.crs)

apiary.coords.sf <- st_intersection(apiary.aea, st_make_valid(maryland.aea))
apiary.coords <- apiary.coords.sf %>%
  st_coordinates()
# Calculate distance between wild bee locations and each apiary
curr.dist <- rdist(data.list$coords, apiary.coords)
# Set distances tr NA if > 3.5
curr.dist <- ifelse(curr.dist > 3.5, NA, curr.dist)
# Total number of apiaries within 3.5km of each wild bee survey location
n.apiary.within.range <- apply(curr.dist, 1, function(a) sum(!is.na(a)))
n.apiary.within.range

# Generate plot showing relationship between the apiary density metric used
# in the manuscript and a static 3.5km density metric. 
plot.df <- data.frame(apiary = data.list$covs$apiary.density, 
                      comparison = n.apiary.within.range)

comparison.metric.plot <- ggplot(data = plot.df, aes(x = comparison, y = apiary)) + 
  geom_point() + 
  labs(x = 'Number of apiaries within 3.5km', y = 'Apiary density metric', title = '(C)') + 
  theme_bw(base_size = 16) + 
  theme(text = element_text(family = 'LM Roman 10'), 
        plot.title = element_text(size = 16),
        panel.grid = element_blank()) 

# Figure S1 ---------------------------------------------------------------
# ggarrange(inv.exp.plot, example.apiary.plot, ncol = 2, widths = c(0.4, 0.6))
(inv.exp.plot + example.apiary.plot) / comparison.metric.plot 
ggsave(file = 'figures/Figure-S1.png', device = 'png', 
       units = 'in', width = 10, height = 10)
  