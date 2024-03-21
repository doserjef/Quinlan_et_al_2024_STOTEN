# 1c-get-apiary-surface.R: this script extracts the metric of apiary density
#                          across a 100x100m grid of Maryland and generates 
#                          Figure 1 in the associated manuscript.
# Author: Jeffrey W. Doser and Gabriela M. Quinlan
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(stars)
library(fields)
library(ggpubr)
library(ggthemes)
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
  labs(x = 'Longitude', y = 'Latitude', title = '(A) Apiary density', 
       fill = 'Density') + 
  theme(legend.position = c(0.3, 0.31), 
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
  labs(x = 'Longitude', y = 'Latitude', title = '(B) Landcover', 
       fill = 'Type') + 
  theme(legend.position = c(0.3, 0.31), 
        legend.background = element_rect(fill = NA), 
        text = element_text(family = 'LM Roman 10'),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))
ggarrange(apiary.plot, cdl.plot, ncol = 1)
ggsave(file = 'figures/Figure-1.png', units = 'in', width = 8, height = 8, bg = 'white')
