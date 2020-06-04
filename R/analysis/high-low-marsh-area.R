# How much low marsh
# How much high marsh
# Uncertainty

library(foreign)
library(tidyverse)

midAtlantic <- as_tibble(read.dbf("data/probLowMarsh-dbfs/probLowMarsh_MidAtlantic.tif.vat.dbf",
                   as.is = T)) %>% 
  mutate(region = "mid-Atlantic")

attach(midAtlantic)
plot(x=Value,
     y=Count, type="l")

northeast <- as_tibble(read.dbf("data/probLowMarsh-dbfs/probLowMarsh_Northeast.tif.vat.dbf",
                                  as.is = T)) %>% 
  mutate(region = "Northeast")

attach(northeast)
plot(x=Value,
     y=Count, type="l")

northwest <- as_tibble(read.dbf("data/probLowMarsh-dbfs/probLowMarsh_Northwest.tif.vat.dbf",
                                as.is = T)) %>% 
  mutate(region = "Northwest")

attach(northwest)
plot(x=Value,
     y=Count, type="l")

southCentral <- as_tibble(read.dbf("data/probLowMarsh-dbfs/probLowMarsh_SouthCentral.tif.vat.dbf",
                                as.is = T)) %>% 
  mutate(region = "South Central")

attach(southCentral)
plot(x=Value,
     y=Count, type="l")

southEast <- as_tibble(read.dbf("data/probLowMarsh-dbfs/probLowMarsh_Southeast.tif.vat.dbf",
                                   as.is = T)) %>% 
  mutate(region = "Southeast")

attach(southEast)
plot(x=Value,
     y=Count, type="l")

southWest <- as_tibble(read.dbf("data/probLowMarsh-dbfs/probLowMarsh_Southwest.tif.vat.dbf",
                                as.is = T)) %>% 
  mutate(region = "Southwest")

attach(southWest)
plot(x=Value,
     y=Count, type="l")

allRegions <- midAtlantic %>% 
  bind_rows(northeast) %>% 
  bind_rows(northwest) %>% 
  bind_rows(southCentral) %>% 
  bind_rows(southEast) %>% 
  bind_rows(southWest) %>% 
  mutate(Value = Value / 100)

regionalSummaries <- allRegions %>% 
  group_by(region) %>% 
  summarise(mean_low_marsh = sum(Value * Count),
            mean_high_marsh = sum((1-Value)*Count),
            standard_error = sqrt(sum(Count*Value*(1-Value)))) %>% 
  mutate(pct_low = (mean_low_marsh / (mean_low_marsh+mean_high_marsh)) * 100,
         mean_low_marsh = mean_low_marsh * 900 * 0.0001,
         mean_high_marsh = mean_high_marsh * 900 * 0.0001,
         standard_error = standard_error  * 900 * 0.0001)

nationalSummaries <- allRegions %>% 
  group_by() %>% 
  summarise(mean_low_marsh = sum(Value * Count),
            mean_high_marsh = sum((1-Value)*Count),
            standard_error = sqrt(sum(Count*Value*(1-Value)))) %>% 
  mutate(pct_low = (mean_low_marsh / (mean_low_marsh+mean_high_marsh)) * 100,
         mean_low_marsh = mean_low_marsh * 900 * 0.0001,
         mean_high_marsh = mean_high_marsh * 900 * 0.0001,
         standard_error = standard_error  * 900 * 0.0001) %>% 
  mutate(region = "CONUS")

totalSummaries <- regionalSummaries %>% 
  bind_rows(nationalSummaries) %>% 
  mutate(region = factor(region, 
                         levels = c("Northwest",
                                    "Southwest",
                                    "South Central",
                                    "Southeast",
                                    "mid-Atlantic",
                                    "Northeast",
                                    "CONUS"))) %>% 
  arrange(region)

write_csv(totalSummaries,
          "tables/high_low_marsh_summaries.csv")
