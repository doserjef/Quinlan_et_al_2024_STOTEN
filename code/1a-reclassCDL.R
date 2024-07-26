# 1a-reclassCDL.R: script to reclass the Cropland Data Layer into the five categories
#                  shown in Figure 1 that were used in subsequent analyses.
#                  NOTE: the raw CDL data are not included on GitHub because of file
#                        size limitations and so this script will not run successfully. 
#                        If the raw CDL files are desired, they can be downloaded from 
#                        cropscape or contact the corresponding author.
# Author: Gabriela M. Quinlan
library(raster); library(foreign); library(dplyr)

# downloaded from cropscape on 10-24-23
# WGS84 projection of Maryland 2018

# Load in rasters for Maryland and DC area.
mdDCclip <- raster('data/CDL_2018_clip_20231027224359_779818224.tif', RAT = TRUE)
mdDCdbf <- read.dbf('data/CDL_2018_clip_20231027224359_779818224.tif.vat.dbf')

# Specify attribute levels in the CDL for reclassing ----------------------
other_crops <- c("Alfalfa",  "Almonds", "Apples", "Apricots", "Asparagus", 
                 "Avocados", "Blueberries", "Broccoli", "Cabbage", 
                 "Caneberries", "Canola", "Cantaloupes", "Carrots", 
                 "Cauliflower", "Celery", "Cherries",  "Chick Peas", 
                 "Christmas Trees", "Citrus",  "Cotton",
                 "Cranberries", "Cucumbers", "Dbl Crop Lettuce/Cantaloupe", 
                 "Dbl Crop Lettuce/Cotton", "Dry Beans", "Eggplants", 
                 "Garlic", "Gourds",  "Grapes", "Greens", "Honeydew Melons", 
                 "Herbs", "Hops", "Lentils", "Lettuce", "Mint", 
                 "Misc Vegs & Fruits" , "Mustard", "Nectarines", 
                 "Olives", "Onions", "Oranges", "Other Crops" , "Other Tree Crops",
                 "Peaches", "Peanuts" , "Pears" ,"Peas", "Pecans",  
                 "Peppers", "Pistachios", "Plums",  "Pomegranates", 
                 "Potatoes",  "Prunes" ,"Pumpkins" , "Radishes",  "Rape Seed", 
                 "Safflower", "Sod/Grass Seed" , "Squash", 
                 "Strawberries", "Sugarbeets", "Sugarcane", "Sunflower", 
                 "Sweet Potatoes" ,  "Tobacco", "Tomatoes" ,  "Turnips",
                 "Walnuts",  "Watermelons")
cornSoyGrain <- c("Barley", "Buckwheat", "Camelina", "Corn", "Dbl Crop Barley/Corn", 
                  "Dbl Crop Barley/Sorghum", "Dbl Crop Barley/Soybeans",
                  "Dbl Crop Corn/Soybeans", "Dbl Crop Durum Wht/Sorghum", 
                  "Dbl Crop Lettuce/Barley", "Dbl Crop Lettuce/Durum Wht",  
                  "Dbl Crop Oats/Corn",  "Dbl Crop Soybeans/Cotton" , 
                  "Dbl Crop Soybeans/Oats" ,  "Dbl Crop Triticale/Corn",  
                  "Dbl Crop WinWht/Corn", "Dbl Crop WinWht/Cotton", 
                  "Dbl Crop WinWht/Sorghum",  "Dbl Crop WinWht/Soybeans",  
                  "Durum Wheat" , "Flaxseed", "Millet", "Oats",  
                  "Other Hay/Non Alfalfa", "Other Small Grains",  "Pop or Orn Corn", 
                  "Rice", "Rye",  "Sorghum", "Soybeans" , "Speltz", 
                  "Spring Wheat" , "Sweet Corn", "Triticale", "Winter Wheat")
pasture <- c("Clover/Wildflowers", "Fallow/Idle Cropland",  
             "Grassland/Pasture", "Herbaceous Wetlands", "Switchgrass", "Vetch")

na <- c("Barren", "Clouds/No Data", "Aquaculture","Background", 
        "Nonag/Undefined","Open Water", "Perennial Ice/Snow" , 
        "Water", "Wetlands")

forest <- c("Deciduous Forest", "Evergreen Forest", "Forest", 
            "Mixed Forest", "Shrubland", "Woody Wetlands")

developed <- c("Developed", "Developed/High Intensity", "Developed/Low Intensity", 
               "Developed/Med Intensity", "Developed/Open Space")

# Key for reclassing to five categories for use in analysis
reclassKey <- as.data.frame(cbind(levels(dbf$CLASS_NAME), rep(NA, length(levels(dbf$CLASS_NAME))))) %>%
  rename("CLASS_NAME" = "V1", 
         "reclass" = "V2") %>%
  mutate(reclass = ifelse(CLASS_NAME %in% c(other_crops, cornSoyGrain), "crop", reclass)) %>%
  #mutate(reclass = ifelse(CLASS_NAME %in% c(cornSoyGrain), "cornSoy", reclass)) %>%
  mutate(reclass = ifelse(CLASS_NAME %in% c(pasture), "pasture", reclass)) %>%
  mutate(reclass = ifelse(CLASS_NAME %in% c(na), "na", reclass)) %>%
  mutate(reclass = ifelse(CLASS_NAME %in% c(forest), "forest", reclass)) %>%
  mutate(reclass = ifelse(CLASS_NAME %in% c(developed), "developed", reclass)) %>%
  full_join(dbf) %>% 
  na.omit() %>%
  mutate(reclassVal = as.numeric(as.factor(reclass))) 

# Reclass the raster based on above delineations --------------------------  
reclassMDDC <- reclassify(mdDCclip, reclassKey[,c("VALUE", "reclassVal")])  

# Save reclassed raster to hard drie --------------------------------------
writeRaster(reclassMDDC,'data/reclassMDDC.tif')
