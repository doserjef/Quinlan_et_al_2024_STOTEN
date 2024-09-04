# 7b-get-apiary-surface-map.R: this script extracts the metric of apiary density
#                              across a 100x100m grid of Maryland and generates 
#                              Figure 1 in the associated manuscript. This script
#                              also generates histograms to showcase the spread of
#                              the observed covariate data relative to the covariate
#                              data across the state. 
# Author: Jeffrey W. Doser and Gabriela M. Quinlan
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(stars)
library(fields)
library(patchwork)
library(RColorBrewer)

# NOTE: the apiary locations are propietary and not available on GitHub, 
#       and thus this script will not run successfully. 

# Get prediction coordinates ----------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
maryland <- usa %>% filter(ID == 'maryland')
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
maryland <- maryland %>%
  st_transform(my.crs)

# Grid the area for predictions. 
# dx and dy are in km as above. Can change as desired.
grid.pred <- st_as_stars(st_bbox(maryland), dx = 1, dy = 1)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
                           coords = c('x', 'y'), 
                           crs = st_crs(maryland))

# Intersect with region of interest
# Takes a few minutes
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(maryland))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))
# Clear things out to save some memory
rm(coords.pred.sf)
gc()
coords.pred.sf <- st_as_sf(coords.0, 
                           coords = c('X', 'Y'),
                           crs = st_crs(maryland))

# Generate plot for apiary data -------------------------------------------
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

apiary.coords.sf <- st_intersection(apiary.aea, st_make_valid(maryland))
apiary.coords <- apiary.coords.sf %>%
  st_coordinates()
# Calculate distance between grid cell and each apiary
curr.dist <- rdist(coords.0, apiary.coords)
# Set distances to NA if > 6
curr.dist <- ifelse(curr.dist > 6, NA, curr.dist)
# Calculate density metric for each apiary within 6km of a given site.
density.by.site <- t(apply(curr.dist, 1, function(a) 1 / exp(a / 2)))
# Sum up to get apiary density metric for each site
apiary.density <- apply(density.by.site, 1, sum, na.rm = TRUE)

coords.pred.sf$apiary.density <- apiary.density

apiary.plot.df <- data.frame(x = coords.0[, 1],
                             y = coords.0[, 2],
                             apiary.density = apiary.density)
apiary.plot.df <- st_as_stars(apiary.plot.df, dims = c('x', 'y'))
apiary.plot <- ggplot() +
  geom_stars(data = apiary.plot.df, aes(fill = apiary.density), interpolate = TRUE) +
  geom_sf(data = maryland, col = 'black', alpha = 0) +
  scale_fill_viridis_c(na.value = NA, labels = c('Low', 'High'), breaks = c(0, 30)) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Longitude', y = 'Latitude', title = '(B) Apiary density', 
       fill = 'Density') + 
  theme(legend.position = c(0.25, 0.31), 
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

# Generate plot for reclassed CDL -----------------------------------------
# Clear out memory
rm(grid.pred, coords.pred, coords.pred.sf, density.by.site, curr.dist)
gc()
grid.pred <- st_as_stars(st_bbox(maryland), dx = 0.25, dy = 0.25)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
                           coords = c('x', 'y'), 
                           crs = st_crs(maryland))
# Intersect with region of interest
# Takes a few minutes
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(maryland))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))
reclass.cdl <- raster("data/reclassMDDC.tif")
# Maryland cdl values
maryland.cdl <- raster::extract(reclass.cdl, coords.pred.sf)
plot.df <- data.frame(x = coords.0[, 1],
                      y = coords.0[, 2],
                      cdl = factor(maryland.cdl, levels = c(1, 2, 3, 4, 5),
                                   labels = c('Agriculture', 'Developed', 
                                              'Forest', 'NA', 'Grassy/Herbaceous')))
# Developed -> black
# Forest -> green
# Pasture -> yellow
# NA -> blue
# Ag -> orange
my.cols <- c('#E69F00', 'black', '#009E73', '#56B4E9', '#F0E442')
plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
cdl.plot <- ggplot() +
  geom_stars(data = plot.df, aes(fill = cdl)) +
  geom_sf(data = maryland, col = 'black', alpha = 0) +
  scale_fill_manual(values = my.cols, na.value = NA, na.translate = FALSE) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Longitude', y = 'Latitude', title = '(C) Landcover', 
       fill = 'Type') + 
  theme(legend.position = c(0.35, 0.31), 
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

# Generate plot of the wild bee data locations ----------------------------
load('data/spAbundance-data.rda')
coords.sf <- st_as_sf(as.data.frame(data.list$coords), 
                      coords = c('X', 'Y'),
                      crs = my.crs)
points.plot <- ggplot() +
  geom_sf(data = maryland, col = 'black', fill = 'white') +
  geom_sf(data = coords.sf, col = 'black') +
  theme_bw(base_size = 18) + 
  labs(x = 'Longitude', y = 'Latitude', title = '(A) Wild bee collection locations') + 
  theme(legend.position.inside = c(0.3, 0.31), 
        legend.position = 'inside',
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

points.plot / apiary.plot / cdl.plot
ggsave(file = 'figures/Figure-1.png', units = 'in', width = 7, height = 12, bg = 'white')

# Get smaller set of prediction coordinates for calculating devel ---------
# Grid the area for predictions. 
# dx and dy are in km as above. Can change as desired.
grid.pred <- st_as_stars(st_bbox(maryland), dx = 5, dy = 5)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
                           coords = c('x', 'y'), 
                           crs = st_crs(maryland))

# Intersect with region of interest
# Takes a few minutes
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(maryland))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))
# Clear things out to save some memory
rm(coords.pred.sf)
gc()
coords.pred.sf <- st_as_sf(coords.0, 
                           coords = c('X', 'Y'),
                           crs = st_crs(maryland))

# Calculate apiary density at the smaller grid points ---------------------
# Calculate distance between grid cell and each apiary
curr.dist <- rdist(coords.0, apiary.coords)
# Set distances to NA if > 6
curr.dist <- ifelse(curr.dist > 6, NA, curr.dist)
# Calculate density metric for each apiary within 6km of a given site.
density.by.site <- t(apply(curr.dist, 1, function(a) 1 / exp(a / 2)))
# Sum up to get apiary density metric for each site
apiary.density <- apply(density.by.site, 1, sum, na.rm = TRUE)

coords.pred.sf$apiary.density <- apiary.density

# Calculate developed land land cover -------------------------------------
# Key:
# 1 -> Ag
# 2 -> developed
# 3 -> forest
# 4 -> NA
# 5 -> pasture

# Takes a while to run, so can load in results at the bottom of this section
# Extract mean values at point locations
# calc.devel <- function(a, na.rm = TRUE) {
#   mean(a %in% c(2), na.rm = na.rm)
# }
# coords.buffered <- st_buffer(coords.pred.sf, dist = 0.5)
# coords.buffered <- coords.buffered %>%
#   st_transform(st_crs(reclass.cdl))
# J <- nrow(coords.buffered)
# prop.devel <- rep(NA, J)
# vals <- split(1:J, ceiling(seq_along(1:J) / 50))
# for (j in 1:length(vals)) {
#   print(paste0("Currently on piece ", j, " out of ", length(vals)))
#   indx <- vals[[j]]
#   prop.devel[indx] <- raster::extract(reclass.cdl, y = coords.buffered[indx, ], fun = calc.devel)
# }
# save(prop.devel, file = 'data/developed-for-violin-plot.rda')
load('data/developed-for-violin-plot.rda')

# Get histograms of data spread -------------------------------------------
my.cols <- viridis(2)
# Apiary
obs.plot.df <- data.frame(apiary = data.list$covs$apiary.density)
pred.plot.df <- data.frame(apiary = apiary.density)
plot.df <- data.frame(apiary = c(data.list$covs$apiary.density, apiary.density), 
                      type = c(rep('Observed Locations', length(data.list$covs$apiary.density)), 
                               rep('Maryland', length(apiary.density))))
apiary.violin.plot <- ggplot(data = plot.df, aes(x = type, y = apiary, fill = type)) + 
  geom_violin(alpha = 0.5) + 
  theme_bw(base_size = 18) + 
  scale_fill_viridis_d() + 
  guides(fill = 'none') + 
  labs(x = 'Location', y = 'Apiary density', title = '(A) Apiary density') + 
  theme(text = element_text(family = 'LM Roman 10'), 
        panel.grid = element_blank())
# Devel
obs.plot.df <- data.frame(devel = data.list$covs$devel)
pred.plot.df <- data.frame(devel = prop.devel)
plot.df <- data.frame(devel = c(data.list$covs$devel, prop.devel), 
                      type = c(rep('Observed Locations', length(data.list$covs$devel)), 
                               rep('Maryland', length(prop.devel))))
devel.violin.plot <- ggplot(data = plot.df, aes(x = type, y = devel, fill = type)) + 
  geom_violin(alpha = 0.5) + 
  theme_bw(base_size = 18) + 
  scale_fill_viridis_d() + 
  guides(fill = 'none') + 
  labs(x = 'Location', y = 'Proportion developed landcover', title = '(B) Developed landcover') + 
  theme(text = element_text(family = 'LM Roman 10'), 
        panel.grid = element_blank())

apiary.violin.plot + devel.violin.plot
ggsave(file = 'figures/Figure-S2.png', units = 'in', width = 12, height = 6, bg = 'white')
