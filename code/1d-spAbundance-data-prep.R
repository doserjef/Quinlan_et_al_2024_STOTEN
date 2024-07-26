# 1d-spAbundance-data-prep.R: takes the "cleanMDOccurEffort.csv" from "1b-mdClean.R" 
#                             and prepares it into format for fitting a JSDM in spAbundance. The 
#                             script also calculates the apiary density metric for 
#                             the wild bee locations in the data set, as well as 
#                             developed land cover.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(maps)
library(sf)
library(raster)
library(stars)
library(fields)

# NOTE: most of this script will not run successfully as many data files are too large 
#       for GitHub or contain propietary information. 

# Format the cleaned data set ---------------------------------------------
occur <- read.csv("data/cleanMDOccurEffort.csv")
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs" 

# Only grab data from Maryland
maryland <- st_as_sf(map('state', region = c('maryland'), 
                         fill = TRUE, plot = FALSE))
maryland.wgs <- maryland %>%
  st_transform(4326)
# Get all unique coordinates from the data
coords.df <- occur %>%
  dplyr::select(x = decimalLongitude, y = decimalLatitude) %>%
  unique()
coords.sf <- st_as_sf(coords.df,
                      coords = c("x", "y"),
                      crs = 4326)
# Some simple plots
# All data in Maryland's bounding box
ggplot() + 
  geom_sf(data = maryland.wgs) + 
  geom_sf(data = coords.sf)
# Only grab coordinates in Maryland
coords.sf.wgs.maryland <- st_filter(coords.sf, maryland.wgs)
ggplot() + 
  geom_sf(data = maryland.wgs) + 
  geom_sf(data = coords.sf.wgs.maryland)
# Reassign Maryland coordinates to coords
coords <- coords.sf.wgs.maryland %>%
  st_coordinates()

# Determine wild bee locations within 6km of a neighboring state ----------
neighbor.states <- st_as_sf(map('state', region = c('pennsylvania', 'district of columbia', 
                                                    'virginia', 'west virginia', 'delaware', 
                                                    'new jersey'), fill = TRUE, plot = FALSE))
neighbor.states.wgs <- neighbor.states %>%
  st_transform(4326)
# Get distances of each point to the boundary of a neighbor state.
vals <- matrix(as.numeric(st_geometry(obj = neighbor.states.wgs) %>% 
                            st_cast(to = 'MULTILINESTRING') %>% 
                            st_distance(y = coords.sf.wgs.maryland)), 
               nrow = nrow(neighbor.states), ncol = nrow(coords.sf.wgs.maryland))
# Filter out sites within 6km of a neighboring state
# Sites within 6km of a neighboring state
bad.sites <- which(apply(vals, 2, function(a) sum(a < 6000)) > 0)
coords <- coords[-bad.sites, ]
coords.sf.wgs.maryland <- coords.sf.wgs.maryland[-bad.sites, ]
coords.sf <- coords.sf[-bad.sites, ]

# Paste columns together for a unique site identifier
coord.id.keep <- paste(coords[, 1], coords[, 2], sep = '-')
occur$coord.id <- paste(occur$decimalLongitude, occur$decimalLatitude, sep = '-')
# Only grab the rows in occur with coordinates in coords
occur <- occur %>%
  filter(coord.id %in% coord.id.keep)

# Extract counts from data ------------------------------------------------
# Wasp families
wasp <- c("", "Agaonidae", "Braconidae", "Chalcididae", "Chrysididae",
          "Crabronidae", "Eumenidae", "Figitidae", "Formicidae", "Ichneumonidae",
          "Larridae", "Leucospidae", "Masaridae", "Mutillidae", "Pompilidae",
          "Pteromalidae", "Scoliidae", "Sphecidae", "Tenthredinidae", "Tiphiidae",
          "Vespidae")
# Extract the observed bee genus
bees <- occur %>%
  filter(order == "Hymenoptera") %>%
  filter(!family %in% wasp) %>%
  filter(genus != "") %>%
  pull(genus) %>%
  unique()
  
# Set the missing genus name to "missing" to ensure columns can be pivoted wider
occur$genus[occur$genus == ''] <- 'missing'
# Format data into counts, filling in 0s where appropriate
# Also filter out pitFallTraps
dat <- occur %>%
  filter(!sampleProtocol %in% c('pitFall', 'panTrap_pitFall')) %>%
  group_by(decimalLatitude, decimalLongitude, startDate, sampleProtocol, 
           day, month, year, genus, nTraps) %>%
  summarize(count = n(), 
            timeEffort = unique(timeEffort)) %>%
  pivot_wider(names_from = genus, values_from = count, values_fill = 0) %>%
  group_by(decimalLatitude, decimalLongitude) %>%
  mutate(site = as.factor(cur_group_id())) %>%
  ungroup() 
# Format for spAbundance --------------------------------------------------
# Format coordinates ------------------
coords.df <- dat %>%
  dplyr::select(x = decimalLongitude, y = decimalLatitude) %>%
  unique()
coords.sf <- st_as_sf(coords.df,
                      coords = c("x", "y"),
                      crs = 4326)
# Convert to aea
coords.aea <- coords.sf %>%
  st_transform(my.crs)
coords <- st_coordinates(coords.aea)
# Number of sites
J <- n_distinct(dat$site)
# Number of replicate surveys at each site
K <- dat %>%
  group_by(site) %>%
  summarize(val = n()) %>%
  pull(val)
K.max <- max(K)
# Index of columns that hold genus data
tmp <- dat %>%
  dplyr::select(-decimalLongitude, -decimalLatitude, -startDate, -sampleProtocol, 
                -day, -month, -year, -site, -timeEffort)
# Genus names
sp.names <- names(tmp)
# Number of unique genus. 
N <- length(sp.names)
# Array to hold detection-nondetection data.
y <- array(NA, dim = c(N, J, K.max))
# Matrix of the day the sampling took place
day.effect <- matrix(NA, J, K.max)
# Matrix of sample protocol 
sample.protocol <- matrix(NA, J, K.max)
# Matrix of time effort. Note the time effort variable is the amount of time spent sampling
# for any given survey. For pan traps, it is the amount of time spent surveying times the 
# total number of traps, which represents an overall metric of "sampling effort".
time.effort <- matrix(NA, J, K.max)
# Matrix of year each survey took place
year <- matrix(NA, J, K.max)
# Fill in data
for (j in 1:J) {
  # One site at a time
  tmp <- dat %>%
    filter(site == j)
  # Get the count data for all genera 
  tmp.sp <- tmp %>%
    dplyr::select(-decimalLongitude, -decimalLatitude, -startDate, -sampleProtocol, 
  	 -day, -month, -year, -site, -timeEffort)
  # Day effect
  day.effect[j, 1:K[j]] <- as.numeric(paste(tmp$day, tmp$month, tmp$year, sep = ''))
  # Year
  year[j, 1:K[j]] <- tmp$year
  # Sample protocol
  sample.protocol[j, 1:K[j]] <- tmp$sampleProtocol
  # Get effort (time for non pan trap, time * nTraps for panTrap)
  time.effort[j, 1:K[j]] <- ifelse(tmp$sampleProtocol == 'panTrap', 
				   tmp$timeEffort * ifelse(is.na(tmp$nTraps), 1, tmp$nTraps),
				   tmp$timeEffort)
  y[, j, 1:K[j]] <- t(as.matrix(tmp.sp))
} # j (site)
rownames(y) <- sp.names

# Standardize the amount of effort within each protocol -------------------
hand.net.indx <- which(sample.protocol == 'handNet', arr.ind = TRUE)
pan.trap.indx <- which(sample.protocol == 'panTrap', arr.ind = TRUE)
vane.trap.indx <- which(sample.protocol == 'vaneTrap', arr.ind = TRUE)
# Set the two incorrect time.effort values for handNet to 0
time.effort[hand.net.indx][which(time.effort[hand.net.indx] == 1440)] <- 0
# Set all 0 values to 1 to indicate a very short amount of sampling
time.effort[which(time.effort == 0)] <- 1
# Get the log of time.effort so the estimated relationship with bee counts is linear
time.effort <- log(time.effort)
# Get means and standard deviations
mean.effort.pan.trap <- mean(time.effort[pan.trap.indx])
sd.effort.pan.trap <- sd(time.effort[pan.trap.indx])
mean.effort.hand.net <- mean(time.effort[hand.net.indx])
sd.effort.hand.net <- sd(time.effort[hand.net.indx])
mean.effort.vane.trap <- mean(time.effort[vane.trap.indx])
sd.effort.vane.trap <- sd(time.effort[vane.trap.indx])
time.effort[pan.trap.indx] <- (time.effort[pan.trap.indx] - mean.effort.pan.trap) / sd.effort.pan.trap
time.effort[hand.net.indx] <- (time.effort[hand.net.indx] - mean.effort.hand.net) / sd.effort.hand.net
time.effort[vane.trap.indx] <- (time.effort[vane.trap.indx] - mean.effort.vane.trap) / sd.effort.vane.trap

# Filter out wasps and missing genus --------------------------------------
good.sp <- which(sp.names %in% bees)
y <- y[good.sp, , ]
sp.names <- dimnames(y)[[1]]
# Look at total number obs individuals per genus
sum.vals <- apply(y, 1, sum, na.rm = TRUE)
# Remove species with less than 10 observations
bad.sp <- which(sum.vals < 10)
y <- y[-bad.sp, , ]
sp.names <- dimnames(y)[[1]]


# Create spAbundance data list --------------------------------------------
# Remove Apis
apis.indx <- which(sp.names == 'Apis')
y <- y[-apis.indx, , ]
sp.names <- dimnames(y)[[1]]
# Get into data list for spAbundance
data.list <- list(y = y, 
                  covs = list(day.random = day.effect, 
                              protocol = sample.protocol, 
                              effort = time.effort, 
                              year = year), 
                  coords = coords)
# Reformat y to improve model fit by putting a common species first
start.sp <- c('Lasioglossum', 'Megachile', 'Halictus', 'Ceratina', 'Melitoma')
# Other species code
indices <- rep(NA, 4)
for (i in 1:length(start.sp)) {
  indices[i] <- which(sp.names == start.sp[i])
}
indices.other <- 1:nrow(data.list$y)
indices.other <- indices.other[-indices]
# Ordered y
y.ordered <- data.list$y[c(indices, indices.other), , ]
data.list$y <- y.ordered

# Get apiary density at the point locations -------------------------------
# Get ccoordinates in AEA
coords.buffered <- st_buffer(coords.aea, dist = 6)
# Check to make sure the complete circle is within Maryland boundaries
maryland.aea <- st_transform(maryland.wgs, my.crs)
# Check to make sure all buffered points fall completely within Maryland
ggplot() + 
  geom_sf(data = maryland.aea) + 
  geom_sf(data = coords.buffered)
# Load apiary data --------------------------------------------------------
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
curr.dist <- rdist(coords, apiary.coords)
# Set distances tr NA if > 6
curr.dist <- ifelse(curr.dist > 6, NA, curr.dist)
# Total number of apiaries within 6km of each wild bee survey location
n.apiary.within.6km <- apply(curr.dist, 1, function(a) sum(!is.na(a)))
n.apiary.within.6km
# Mean distance to the apiaries within 6km
mean.dist.to.apiary <- apply(curr.dist, 1, mean, na.rm = TRUE)
mean.dist.to.apiary
# Calculate density metric for each apiary within 6km of a given site.
density.by.site <- t(apply(curr.dist, 1, function(a) 1 / exp(a / 2)))
# Sum up to get apiary density metric for each site
apiary.density <- apply(density.by.site, 1, sum, na.rm = TRUE)

coords.aea$apiary.density <- apiary.density
# Simple plot of apiary density at points
ggplot() + 
  geom_sf(data = maryland.aea) + 
  geom_sf(data = coords.aea, aes(col = apiary.density), size = 3) +
  scale_color_viridis_c()

# Plot of where the actual apiaries are, and can explore the number of 
# apiaries within a given boundary of a wild bee location
ggplot() + 
  geom_sf(data = maryland.aea) + 
  geom_sf(data = coords.buffered[62, ], fill = 'white') + 
  geom_sf(data = apiary.coords.sf, col = 'red')

data.list$covs$apiary.density <- apiary.density 

# Calculate four landcover variables around apiary locations --------------
# Key:
# 1 -> Ag
# 2 -> developed
# 3 -> forest
# 4 -> NA
# 5 -> pasture

# Load in reclassed CDL
reclass.cdl <- raster("data/reclassMDDC.tif")

# Extract mean values at point locations
calc.ag <- function(a, na.rm = TRUE) {
  mean(a %in% c(1), na.rm = na.rm)
}
calc.devel <- function(a, na.rm = TRUE) {
  mean(a %in% c(2), na.rm = na.rm)
}
calc.forest <- function(a, na.rm = TRUE) {
  mean(a %in% c(3), na.rm = na.rm)
}
calc.pasture <- function(a, na.rm = TRUE) {
  mean(a %in% c(5), na.rm = na.rm)
}
calc.NA <- function(a, na.rm = TRUE) {
  mean(a %in% c(4), na.rm = na.rm)
}
my.length <- function(a, na.rm = FALSE) {
  sum(!is.na(a))
}
coords.aea.buffered <- st_buffer(coords.aea, dist = 0.5)
coords.aea.buffered <- coords.aea.buffered %>%
  st_transform(st_crs(reclass.cdl))
J <- nrow(coords.aea)
prop.ag <- rep(NA, J)
prop.devel <- rep(NA, J)
prop.forest <- rep(NA, J)
prop.pasture <- rep(NA, J)
prop.NA <- rep(NA, J)
vals <- split(1:J, ceiling(seq_along(1:J) / 10))
for (j in 1:length(vals)) {
  print(paste0("Currently on piece ", j, " out of ", length(vals)))
  indx <- vals[[j]]
  prop.NA[indx] <- raster::extract(reclass.cdl, y = coords.aea.buffered[indx, ], fun = calc.NA)
  prop.ag[indx] <- raster::extract(reclass.cdl, y = coords.aea.buffered[indx, ], fun = calc.ag)
  prop.devel[indx] <- raster::extract(reclass.cdl, y = coords.aea.buffered[indx, ], fun = calc.devel)
  prop.forest[indx] <- raster::extract(reclass.cdl, y = coords.aea.buffered[indx, ], fun = calc.forest)
  prop.pasture[indx] <- raster::extract(reclass.cdl, y = coords.aea.buffered[indx, ], fun = calc.pasture)
}

covs.prop.land <- data.frame(ag = prop.ag,
                             devel = prop.devel,
                             forest = prop.forest,
                             other = prop.NA,
                             pasture = prop.pasture)
covs.df <- data.frame(covs.prop.land,
                      apiary.density = apiary.density)
# Devel and apiary are fairly positively correlated
cor(covs.df)

data.list$covs$devel <- prop.devel
# Save results in data list for spAbundance -------------------------------
save(data.list, file = 'data/spAbundance-data.rda')

# Save data in flat file if desired
dat.save <- dat %>%
  dplyr::select(decimalLatitude, decimalLongitude, startDate, sampleProtocol, 
                day, month, year, timeEffort, site, any_of(dimnames(data.list$y)[[1]]))
write.csv(dat.save, file = 'data/data-spAbundance-flat.csv', row.names = FALSE)
