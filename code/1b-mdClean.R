# 1b-mdClean.R: script to clean the Marlyand data from Sam Droege and colleagues 
#            at the USGS. 
# Author: Gabriela M. Quinlan and Melanie Kammerer
rm(list = ls())
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

# NOTE: these raw data are too large to include on GitHub, and thus this script will not
#       successfully run. These data can be provided upon request to the corresponding 
#       author, or can be downloaded directly from 
#       https://data.usgs.gov/datacatalog/data/USGS:29e35741-5920-4563-97fb-d555a9ca24ed
occur <- read.table("gbif-download-oct-21-2023/occurrence.txt", sep='\t', header=T, quote='', comment='')
# downloaded 10-21-23

# excluded protocols: "" (unknown); "malaise trap" (n=15)

# Filter for pan trap data ------------------------------------------------
panTrap <- occur %>% 
  select_if(~sum(!is.na(.)) > 0) %>%
  # Restrict to 2017 and 2018
  filter(year == 2017 | year== 2018) %>%
  # Restrict to Maryland
  filter(decimalLatitude > 37.8713 & decimalLatitude < 39.7425) %>% # maryland bounding box
  filter(decimalLongitude > -79.4938 & decimalLongitude < -75.0450) %>%
  # Only pan traps 
  filter(samplingProtocol == "pan trap") %>%
  # Only grab columns of interest
  select(eventDate, year, month, day, verbatimEventDate, fieldNotes, sampleSizeValue, eventRemarks, decimalLatitude, decimalLongitude, stateProvince, kingdom, phylum, class, order, family, genus) %>% 
  # Get columns to aid in sampling effort determination
  tidyr::separate(col=eventRemarks, sep=";", into=c('TrapVolume','TrapColor','TrapLiquid'))  %>%
  # Get time spent sampling and store in timeEffort column
  mutate(verbatimEventDate = str_replace_all(verbatimEventDate, pattern = "[x]", replacement = "0")) %>%
  tidyr::separate(col=verbatimEventDate, into=c("time1", 'time2'), sep=";") %>%
  tidyr::separate(col=time1, sep="StartDateTime:", into=c("delete", "time1")) %>%
  mutate(startDate = as.POSIXct(substring(time1, 1, 12), format = "%Y%m%d%H%M")) %>%  # remainder are times and fillers
  tidyr::separate(col=time2, sep="EndDateTime:", into=c("delete", "time2")) %>%
  mutate(endDate = as.POSIXct(substring(time2, 1, 12), format = "%Y%m%d%H%M")) %>% # remainder are times and fillers 
  mutate(timeEffort = endDate - startDate) %>%
  # Fix samples where the year is incorrectly recorded.
  mutate(timeEffort = ifelse(timeEffort >= 31536000, timeEffort - 31536000, timeEffort)) %>% # 1 yr, likely year is wrong
  mutate(timeEffort = timeEffort / 60) %>% # convert to minutes
  # Determine ratios of blue, white, and yellow traps
  mutate(blue = NA, 
         white = NA, 
         yellow = NA) %>%
  mutate(blue = ifelse(grepl("10 blue", fieldNotes), 10, blue), 
         white = ifelse(grepl("10 white", fieldNotes), 10, white),
         yellow = ifelse(grepl("10 yellow", fieldNotes), 10, yellow)) %>%
  mutate(blue = ifelse(grepl("18 blue", fieldNotes), 18, blue), 
         white = ifelse(grepl("18 white", fieldNotes), 18, white),
         yellow = ifelse(grepl("18 yellow", fieldNotes), 18, yellow)) %>%
  mutate(blue = ifelse(grepl("36 blue", fieldNotes), 36, blue), 
         white = ifelse(grepl("36 white", fieldNotes), 36, white),
         yellow = ifelse(grepl("36 yellow", fieldNotes), 36, yellow)) %>%
  mutate(blue = ifelse(grepl("8 blue", fieldNotes), 8, blue), 
         white = ifelse(grepl("8 white", fieldNotes), 8, white),
         yellow = ifelse(grepl("8 yellow", fieldNotes), 8, yellow)) %>%
  mutate(blue = ifelse(grepl("8 fl b", fieldNotes), 8, blue), 
         white = ifelse(grepl("8 w", fieldNotes), 8, white),
         yellow = ifelse(grepl("8 fl y", fieldNotes), 8, yellow)) %>%
  mutate(yellow = ifelse(grepl("8 fl yl", fieldNotes), 8, yellow)) %>%
  mutate(blue = ifelse(grepl("9FBL", fieldNotes), 9, blue), 
         white = ifelse(grepl("9W", fieldNotes), 9, white),
         yellow = ifelse(grepl("9FYL", fieldNotes), 9, yellow)) %>%
  mutate(blue = ifelse(grepl("9FLB", fieldNotes), 9, blue)) %>%
  mutate(blue = ifelse(grepl("3 blue", fieldNotes), 3, blue), 
         white = ifelse(grepl("3 white", fieldNotes), 3, white),
         yellow = ifelse(grepl("3 yellow", fieldNotes), 3, yellow)) %>%
  mutate(blue = ifelse(grepl("10 fluorescent blue", fieldNotes), 10, blue), 
         white = ifelse(grepl("10 fluorescent blue and white bowls", fieldNotes), 10, white),
         yellow = ifelse(grepl("10 fluorescent yellow", fieldNotes), 10, yellow)) %>%
  mutate(blue = ifelse(grepl("10flbl", fieldNotes), 10, blue), 
         white = ifelse(grepl("10white", fieldNotes), 10, white),
         yellow = ifelse(grepl("10 flyl", fieldNotes), 10, yellow)) %>%
  mutate(yellow = ifelse(grepl("9FLYL", fieldNotes), 9, yellow)) %>%
  # Determine trap volumes
  mutate(TrapVolume = ifelse(TrapVolume == "TrapVolume:in field note" & grepl("12 oz beer cups", fieldNotes), 12, TrapVolume)) %>%
  mutate(TrapVolume = ifelse(grepl("3.25oz", TrapVolume), 3.25, TrapVolume)) %>%
  mutate(TrapVolume = ifelse(grepl("12.0oz", TrapVolume), 12, TrapVolume)) %>%
  mutate(TrapVolume = ifelse(grepl("2.0oz", TrapVolume), 2, TrapVolume)) %>% 
  mutate(TrapVolume = ifelse(TrapVolume == "TrapVolume:in field note" & grepl("3.25", fieldNotes), 12, TrapVolume)) %>%
  mutate(TrapLiquid = ifelse(TrapColor == " TrapLiquid:glycol propylene",  " TrapLiquid:glycol propylene", TrapLiquid)) %>%
  mutate(white = ifelse(grepl("Fail: 1 w", fieldNotes), white - 1, white)) %>%
  mutate(white = ifelse(grepl("Fail: 2 q", fieldNotes), white - 2, white)) %>%
  mutate(white = ifelse(grepl("fail: 1 w", fieldNotes), white - 1, white)) %>%
  mutate(white = ifelse(grepl("Fail: 1 b 1 w", fieldNotes), white - 1, white)) %>%
  mutate(blue = ifelse(grepl("Fail: 1 b", fieldNotes), blue - 1, blue)) %>%
  mutate(blue = ifelse(grepl("Fail: 1 w 1 b", fieldNotes), blue - 1, blue)) %>%
  mutate(blue = ifelse(grepl("Fail: 1 b 1 w", fieldNotes), blue - 1, blue)) %>%
  mutate(yellow = ifelse(grepl("Fail: 1 w 1 b 1 y", fieldNotes), yellow - 1, yellow)) %>%
  mutate(yellow = ifelse(grepl("Fail: 1 y", fieldNotes), yellow - 1, yellow)) %>%
  mutate(yellow = ifelse(grepl("I failed y", fieldNotes), yellow - 1, yellow)) %>%
  mutate(nTraps = blue + white + yellow) %>%
  # Manually manipulate unusual instances of traps to get accurate number of functional
  # traps for each observation
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 21", fieldNotes), 21,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 28", fieldNotes), 28,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 29", fieldNotes), 29,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 30", fieldNotes), 30,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 32", fieldNotes), 32,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 33", fieldNotes), 33,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("FTE: Fresh Tidal Estuary;Bee Bowls: 38", fieldNotes), 38,  nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 4 bowls", fieldNotes), nTraps - 4, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("lost 2 bowls", fieldNotes), nTraps - 2, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 6 bowls", fieldNotes), nTraps - 6, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("lost 24 bowls", fieldNotes), nTraps - 24, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 7 bowls", fieldNotes), nTraps - 7, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 3 bowls", fieldNotes), nTraps - 3, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("lost 10 bowls", fieldNotes), nTraps - 10, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 8 bowls", fieldNotes), nTraps - 8, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 19 bowls", fieldNotes), nTraps - 19, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 11 bowls", fieldNotes), nTraps - 11, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Lost 2 bowls", fieldNotes), nTraps - 2, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("1 missing", fieldNotes), nTraps - 1, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("destroying all but 2 of 27 cups", fieldNotes), nTraps - 25, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("1 bowl missing", fieldNotes), nTraps - 1, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("5 bowls flattened so no contents, 15 bowls missing", fieldNotes), nTraps - 20, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("1 bowl crushed", fieldNotes), nTraps - 1, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("1 lost", fieldNotes), nTraps - 1, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("17 bowls shredded", fieldNotes), nTraps - 17, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("1 bowl down", fieldNotes), nTraps - 1, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("9 bowls crushed", fieldNotes), nTraps - 9, nTraps)) %>%
  mutate(nTraps = ifelse(!is.na(sampleSizeValue), sampleSizeValue, nTraps)) %>% # feel like this should trump the other notes 
  dplyr::select(- c(delete, TrapColor)) %>%
  mutate(TrapVolume = as.numeric(TrapVolume)) %>%
  dplyr::select(year, month, day, startDate, endDate, timeEffort, 
                decimalLatitude, decimalLongitude, stateProvince,
                TrapVolume, TrapLiquid, blue, white, yellow, nTraps, 
                class, order, family, genus) %>%
  rename(trapVolOz = TrapVolume, 
         trapLiquid = TrapLiquid) %>%
  mutate(sampleProtocol = "panTrap") # some compared to "50 ml falcon tubes" but unsure if includes these in the count. 
# these also  "only blue and yellow painted ones were used" but did not note any numbers. -- not many though. 
# also did not differentiate "flourecent" from regular blue/ yellow

# Filter handNet data -----------------------------------------------------
handNet <- occur %>% 
  select_if(~sum(!is.na(.)) > 0) %>%
  filter(year == 2017 | year== 2018) %>%
  filter(decimalLatitude > 37.8713 & decimalLatitude < 39.7425) %>% # maryland bounding box
  filter(decimalLongitude > -79.4938 & decimalLongitude < -75.0450) %>%
  filter(samplingProtocol == "hand net") %>%
  select(eventDate, year, month, day, verbatimEventDate, fieldNotes, sampleSizeValue, eventRemarks, decimalLatitude, decimalLongitude, stateProvince, kingdom, phylum, class, order, family, genus, samplingEffort) %>% 
  mutate(verbatimEventDate = str_replace_all(verbatimEventDate, pattern = "[x]", replacement = "0")) %>%
  tidyr::separate(col=verbatimEventDate, into=c("time1", 'time2'), sep=";") %>%
  tidyr::separate(col=time1, sep="StartDateTime:", into=c("delete", "time1")) %>%
  mutate(startDate = as.POSIXct(substring(time1, 1, 12), format = "%Y%m%d%H%M")) %>%  # remainder are times and fillers
  tidyr::separate(col=time2, sep="EndDateTime:", into=c("delete", "time2")) %>%
  mutate(endDate = as.POSIXct(substring(time2, 1, 12), format = "%Y%m%d%H%M")) %>% # remainder are times and fillers 
  mutate(timeEffort = endDate - startDate) %>%
  mutate(timeEffort = ifelse(timeEffort >= 31536000, timeEffort - 31536000, timeEffort)) %>% # 1 yr, likely year is wrong
  mutate(timeEffort = timeEffort / 60) %>% # convert to minutes
  dplyr::select(year, month, day, startDate, endDate, timeEffort, 
                decimalLatitude, decimalLongitude, stateProvince, 
                class, order, family, genus) %>%
  mutate(sampleProtocol = "handNet")

# Filter more pan trap data that are coded up differently -----------------
panTrap2 <- occur %>% 
  select_if(~sum(!is.na(.)) > 0) %>%
  filter(year == 2017 | year== 2018) %>%
  filter(decimalLatitude > 37.8713 & decimalLatitude < 39.7425) %>% # maryland bounding box
  filter(decimalLongitude > -79.4938 & decimalLongitude < -75.0450) %>%
  filter(samplingProtocol == "in field note" & grepl("trap", fieldNotes)) %>%
  select(eventDate, year, month, day, verbatimEventDate, fieldNotes, sampleSizeValue, eventRemarks, decimalLatitude, decimalLongitude, stateProvince, kingdom, phylum, class, order, family, genus) %>% 
  tidyr::separate(col=eventRemarks, sep=";", into=c('TrapVolume','TrapColor','TrapLiquid'))  %>%
  mutate(verbatimEventDate = str_replace_all(verbatimEventDate, pattern = "[x]", replacement = "0")) %>%
  tidyr::separate(col=verbatimEventDate, into=c("time1", 'time2'), sep=";") %>%
  tidyr::separate(col=time1, sep="StartDateTime:", into=c("delete", "time1")) %>%
  mutate(startDate = as.POSIXct(substring(time1, 1, 12), format = "%Y%m%d%H%M")) %>%  # remainder are times and fillers
  tidyr::separate(col=time2, sep="EndDateTime:", into=c("delete", "time2")) %>%
  mutate(endDate = as.POSIXct(substring(time2, 1, 12), format = "%Y%m%d%H%M")) %>% # remainder are times and fillers 
  mutate(timeEffort = endDate - startDate) %>%
  mutate(timeEffort = 60 * as.numeric(timeEffort)) %>% # convert to minutes
  mutate(nTraps = ifelse(grepl("30 3.25 oz bowl traps, 30 100 mL centrifuge tube traps", fieldNotes), 60, NA)) %>%
  mutate(TrapVolume = ifelse(grepl("3.25oz", TrapVolume), 3.25, TrapVolume)) %>%
  mutate(TrapVolume = ifelse(grepl("3.25 oz", fieldNotes), 3.25, TrapVolume)) %>% # 100 ml also ~3.25 oz
  mutate(TrapVolume = ifelse(grepl("solo cup", fieldNotes), 16, TrapVolume)) %>%
  mutate(TrapLiquid = ifelse(grepl("Glycol", fieldNotes), " TrapLiquid:glycol propylene", TrapVolume)) %>%
  mutate(TrapVolume = as.numeric(TrapVolume)) %>%
  mutate(sampleProtocol = ifelse(grepl("30 3.25 oz bowl traps, 30 100 mL centrifuge tube traps", fieldNotes), "panTrap_pitFall", "panTrap")) %>%
  mutate(nTraps = ifelse(!is.na(sampleSizeValue), sampleSizeValue, nTraps)) %>% 
  dplyr::select(year, month, day, startDate, endDate, timeEffort, 
                decimalLatitude, decimalLongitude, stateProvince, 
                TrapVolume, TrapLiquid, nTraps, 
                class, order, family, genus, sampleProtocol) %>%
  rename(trapVolOz = TrapVolume, 
         trapLiquid = TrapLiquid)

# Filter pit fall data ----------------------------------------------------
pitFall <- occur %>% 
  select_if(~sum(!is.na(.)) > 0) %>%
  filter(year == 2017 | year== 2018) %>%
  filter(decimalLatitude > 37.8713 & decimalLatitude < 39.7425) %>% # maryland bounding box
  filter(decimalLongitude > -79.4938 & decimalLongitude < -75.0450) %>%
  filter(samplingProtocol == "pitfall trap") %>% # vane trap
  select(eventDate, year, month, day, verbatimEventDate, fieldNotes, sampleSizeValue, eventRemarks, decimalLatitude, decimalLongitude, stateProvince, kingdom, phylum, class, order, family, genus) %>% 
  tidyr::separate(col=eventRemarks, sep=";", into=c('TrapVolume','TrapColor','TrapLiquid'))  %>%
  mutate(verbatimEventDate = str_replace_all(verbatimEventDate, pattern = "[x]", replacement = "0")) %>%
  tidyr::separate(col=verbatimEventDate, into=c("time1", 'time2'), sep=";") %>%
  tidyr::separate(col=time1, sep="StartDateTime:", into=c("delete", "time1")) %>%
  mutate(startDate = as.POSIXct(substring(time1, 1, 12), format = "%Y%m%d%H%M")) %>%  # remainder are times and fillers
  tidyr::separate(col=time2, sep="EndDateTime:", into=c("delete", "time2")) %>%
  mutate(endDate = as.POSIXct(substring(time2, 1, 12), format = "%Y%m%d%H%M")) %>% # remainder are times and fillers 
  mutate(timeEffort = endDate - startDate) %>%
  mutate(timeEffort = 60 * as.numeric(timeEffort)) %>% # convert to minutes
  mutate(TrapVolume = ifelse(grepl("3.25oz", TrapVolume), 3.25, NA)) %>%
  mutate(nTraps = ifelse(grepl("Bee Bowls: 30", fieldNotes), 30, NA)) %>%
  mutate(nTraps = ifelse(grepl("Bee Bowls: 27", fieldNotes), 27, nTraps)) %>%
  mutate(nTraps = ifelse(grepl("Bee Bowls: 33", fieldNotes), 33, nTraps)) %>%
  mutate(nTraps = ifelse(!is.na(sampleSizeValue), sampleSizeValue, nTraps)) %>% 
  dplyr::select(year, month, day, startDate, endDate, timeEffort, 
                decimalLatitude, decimalLongitude, stateProvince, 
                TrapVolume, TrapLiquid, nTraps, 
                class, order, family, genus) %>%
  rename(trapVolOz = TrapVolume, 
         trapLiquid = TrapLiquid) %>%
  mutate(sampleProtocol = "pitFall")
  
# Filter vane trap data ---------------------------------------------------
vaneTrap <- occur %>% 
  select_if(~sum(!is.na(.)) > 0) %>%
  filter(year == 2017 | year== 2018) %>%
  filter(decimalLatitude > 37.8713 & decimalLatitude < 39.7425) %>% # maryland bounding box
  filter(decimalLongitude > -79.4938 & decimalLongitude < -75.0450) %>%
  filter(samplingProtocol == "vane trap") %>%
  select(eventDate, year, month, day, verbatimEventDate, fieldNotes, sampleSizeValue, eventRemarks, decimalLatitude, decimalLongitude, stateProvince, kingdom, phylum, class, order, family, genus) %>% 
  tidyr::separate(col=eventRemarks, sep=";", into=c('TrapVolume','TrapColor','TrapLiquid'))  %>%
  mutate(verbatimEventDate = str_replace_all(verbatimEventDate, pattern = "[x]", replacement = "0")) %>%
  tidyr::separate(col=verbatimEventDate, into=c("time1", 'time2'), sep=";") %>%
  tidyr::separate(col=time1, sep="StartDateTime:", into=c("delete", "time1")) %>%
  mutate(startDate = as.POSIXct(substring(time1, 1, 12), format = "%Y%m%d%H%M")) %>%  # remainder are times and fillers
  tidyr::separate(col=time2, sep="EndDateTime:", into=c("delete", "time2")) %>%
  mutate(time2 = ifelse(time2 == "20180108160000", "20190108160000", time2)) %>% # negative, looks like wrong year
  mutate(endDate = as.POSIXct(substring(time2, 1, 12), format = "%Y%m%d%H%M")) %>% # remainder are times and fillers 
  mutate(timeEffort = endDate - startDate) %>%
  mutate(timeEffort = as.numeric(timeEffort) * 24 * 60) %>% # convert to minutes
  mutate(TrapColor = ifelse(grepl("blue", TrapColor), "blue", NA)) %>%
  mutate(TrapColor = ifelse(grepl("Blue-Yellow", fieldNotes), "blueYellow", TrapColor)) %>%
  mutate(TrapColor = ifelse(grepl("Yellow Base - Blue Top", fieldNotes), "blueYellow", TrapColor)) %>%
  mutate(TrapColor = ifelse(grepl("Yellow Base-Blue Top", fieldNotes), "blueYellow", TrapColor)) %>%
  mutate(nTraps = ifelse(!is.na(sampleSizeValue), sampleSizeValue, NA)) %>% 
  dplyr::select(year, month, day, startDate, endDate, timeEffort, 
                decimalLatitude, decimalLongitude, stateProvince, 
                TrapLiquid, TrapColor, 
                class, order, family, genus) %>%
  rename(trapColor = TrapColor, 
         trapLiquid = TrapLiquid) %>%
  mutate(sampleProtocol = "vaneTrap")
  
# Combine all sampling approaches in a single data set.
dat <- full_join(panTrap, panTrap2 )%>% 
  full_join(pitFall) %>%
  full_join(vaneTrap) %>%
  full_join(handNet)

# Save cleaned data to hard drive -----------------------------------------
write.csv(dat, "cleanMDOccurEffort.csv", row.names = F)
