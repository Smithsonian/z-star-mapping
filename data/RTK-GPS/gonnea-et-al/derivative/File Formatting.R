# formatting Meagan Gonnea's RTK-GPS data
library(tidyverse)
library(readxl)
library(dplyr)
library(sp)

#devtools::install_github("dkahle/ggmap") # As of time of writing, this is a package that receives frequent 
library(ggmap)
library(rgdal)

# Here's the function itself
convert_UTM_to_latlong <- function(easting, northing, zone) {
  
  # Remove non-numeric characters from zone attribute
  zone <- gsub("[^0-9]", "", zone)
  
  # Change value of digits to ensure as.numeric does not round
  options(digits=22)
  
  # Ensure that easting, northing, and zone are numeric
  easting <- as.numeric(easting)
  northing <- as.numeric(northing)
  zone <- as.numeric(zone)
  
  # Combine the three attributes
  spatialData <- as.matrix(cbind(easting, northing, zone))
  
  # We'll need spatialData as a tibble or data frame going forward
  spatialData <- as_tibble(cbind(easting, northing, zone))
  
  # Now that we flagged the NA rows, we can get rid of them
  spatialData <- na.omit(spatialData)
  
  # And get rid of any zone values that were NA
  zone_list <- as.list(unique(na.omit(zone)))
  
  # Establish the projection we'll convert to. NOTE: this could be a user input
  wgs84 = "+init=epsg:4326"
  
  # Initialize our output dataset
  output <- matrix(nrow = 0, ncol = 2)
  colnames(output) <- c("core_longitude", "core_latitude")
  
  # We'll need to transform the projection separately for data that are in
  #   different zones. So we'll need to subset by the zone values, which we
  #   stored in zone_list
  for (i in 1:length(zone_list)) {
    
    # Filter to just one zone
    spatialData_sub <- spatialData %>%
      filter(zone == zone_list[[i]]) %>%
      select(easting, northing) %>% # We don't need the zone attribute anymore
      na.omit() # Just in case
    
    # Create a dataframe for out subsetted data
    output_sub <- spatialData_sub
    
    # Define the proj4string, using the zone that the subsetted data are in
    proj4string <- CRS(as.character(paste0("+proj=utm +zone=", zone_list[[i]],
                                           " +ellps=WGS84 +datum=WGS84 +units=m +no_defs")))
    
    # Finally, perform the transformation
    sp <-  sp::spTransform(sp::SpatialPoints(list(spatialData_sub$easting, 
                                                  spatialData_sub$northing), proj4string=proj4string),sp::CRS(wgs84))
    
    
    output_sub <- na.omit(sp@coords) # Get rid of NAs again
    colnames(output_sub) <- c("longitude", "latitude") # Rename the output columns
    output <- rbind(output, output_sub) # And slap the data subset into our final output
  }
  
  output <- data.frame(output) # Turn our output into a dataframe
  
  # And output the output
  output
}

# Lat, Long in DD
# m in UTM
# elevation m in NAVD88
# Site
# Year Month Day
# notes, comments, origin ? etc.

# Read in three excel files
sandyNeck <- read_xlsx("data/RTK-GPS/gonnea-et-al/original/AllSandyNeckPoints_RTK2018.xlsx",
                       sheet = 2)
sageLot <- read_xlsx("data/RTK-GPS/gonnea-et-al/original/Sage Lot Pond 2017_2018 RTK Surveys.xlsx",
                     sheet = 1)
restoredMarsh <- read_xlsx("data/RTK-GPS/gonnea-et-al/original/RestoredMarshRTK.xlsx",
                           sheet = 2)
head(restoredMarsh)

restoredMarshEdited <- restoredMarsh %>%
  rename(latitude_dd = lat,
         longitude_dd = long,
         elevation = elev,
         site_id = loc,
         point_id = id,
         transect_id = X__1
         ) %>%
  select(site_id, point_id, latitude_dd, longitude_dd, elevation)

restoredMarshEdited$latitude_dd <- as.numeric(restoredMarshEdited$latitude_dd)

head(sageLot)

sageLotEdited <- sageLot %>%
  mutate(site_id = "Sage Lot") %>%
  rename(elevation = altitude_m)

sageLotEdited$point_id <- as.character(1:nrow(sageLotEdited))

sageLotLatLongs <- convert_UTM_to_latlong(northing = sageLotEdited$northing_m,
                       easting = sageLotEdited$easting_m, zone = "19")

sageLotEdited["latitude_dd"] <- sageLotLatLongs$core_latitude
sageLotEdited["longitude_dd"] <- sageLotLatLongs$core_longitude

sageLotEdited <- sageLotEdited %>%
  mutate(point_id = ifelse(nchar(as.character(point_id)) == 1, paste("SageLot-00", point_id, sep=""),
                                 ifelse(nchar(as.character(point_id)) == 2, paste("SageLot-0", point_id, sep=""), 
                           paste("SageLot-", point_id, sep="")))) %>%
  select(site_id, point_id, latitude_dd, longitude_dd, elevation)


head(sandyNeck)

sandyNeckEdited <- sandyNeck %>%
  rename(easting_m = Easting,
         northing_m = Northing,
         elevation = Elev_NAVD88
         ) %>%
  mutate(site_id = "sandy neck",
         point_id = paste(Descriptor, Point, sep="-"))
  
sandyNeckLotLatLongs <- convert_UTM_to_latlong(northing = sandyNeckEdited$northing_m,
                                            easting = sandyNeckEdited$easting_m, zone = "19" 
  )

sandyNeckEdited["latitude_dd"] <- sandyNeckLotLatLongs$core_latitude
sandyNeckEdited["longitude_dd"] <- sandyNeckLotLatLongs$core_longitude

sandyNeckEdited <- sandyNeckEdited %>% 
  select(site_id, point_id, latitude_dd, longitude_dd, elevation)

sitesRtkCompiled <- bind_rows(sandyNeckEdited, sageLotEdited) %>%
  bind_rows(restoredMarshEdited) %>%
  mutate(site_id = ifelse(site_id == "SGF", "State Game Farm", site_id),
         site_id = ifelse(site_id == "Scusset", "Scussett", site_id)
         ) %>%
  mutate(site_id = toupper(site_id)) %>%
  filter(complete.cases(.))

write_csv(sitesRtkCompiled, "data/RTK-GPS/gonnea-et-al/derivative/GonneaEtAlRtkGps_190702.csv")


MAP_map <- get_stamenmap(bbox = c(left = min(sitesRtkCompiled$longitude_dd, na.rm = T)-0.1,
                                  bottom = min(sitesRtkCompiled$latitude_dd, na.rm = T)-0.1, 
                                  right = max(sitesRtkCompiled$longitude_dd, na.rm = T)+0.1,
                                  top = max(sitesRtkCompiled$latitude_dd, na.rm = T)+0.1))

# Add the CONUS_map with axes set as Longitude and Latitude
ggmap(MAP_map,
      extent = "device",
      xlab = "Longitude",
      ylab = "Latitude") +
  geom_point(data = sitesRtkCompiled, aes(x = longitude_dd,
                                          y = latitude_dd,
                                          color = elevation,
                                          shape = site_id,
                                          size = elevation)) +
  xlab("Longitude") + # Change x axis title
  ylab("Latitude") + # Change y axis title
  scale_shape_manual(values = 0:12) +
  scale_color_gradientn(colors = c("#EFEC35", "#E82323"), na.value = "#1a1a1a")

ggplot(data = sitesRtkCompiled, aes(x = elevation)) +
  geom_histogram() +
  geom_rug(aes(color = site_id), alpha=0.9) +
  scale_color_brewer(palette = "Paired")
