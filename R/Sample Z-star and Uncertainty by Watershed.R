# Need to iterate throught the clipped up versions of Z* and Z* Uncertainty
# Import as a raster, then as a data.frame
# Get summary statistics 
# Median, Mean, s.d., n, Q97.5, 75, 25, 2.5%, min and max
# Add all points to a data frame if under a certain number, randomly sample a subset if above.

library(raster)
library(tidyverse)

# Folder with geospatial files

geo_files <- "F:/z-star-spatial-data/temp/z-star/"

# Get list of geotifs in the specified folder where the geospatial files live
file_list <- list.files(geo_files)
tif_list <- file_list[grepl(".tif", file_list) & !grepl(".tif.aux.xml", file_list) &  !grepl(".tif.xml", file_list)]

maxPoints <- 1000


pb <- txtProgressBar(min = 0, max = length(tif_list), style = 3)

for (i in 1:length(tif_list)) {
  
  input_tif <- tif_list[i]
  
  print(input_tif)
  
  input_path <- paste(geo_files, input_tif, sep="")
  input_raster <- raster(input_path)
  input_df <- as.data.frame(rasterToPoints(input_raster)) 
  
  names(input_df) <- c("x", "y", "z")
  
  input_df <- input_df %>% mutate(Abbrev = str_remove(input_tif, ".tif"))
  
  sum_stats <- input_df %>%
    group_by(Abbrev) %>%
    summarise(mean = mean(z),
              n = n(),
              sd = sd(z),
              min = min(z),
              Q025 = quantile(z, 0.025),
              Q25 = quantile(z,0.25),
              median = median(z),
              Q75 = quantile(z, 0.75),
              Q975 = quantile(z, 0.975),
              max = max(z))
  
  #outlier_thresholds <- sum_stats %>%
  #  mutate(lower_thresh = Q25 - 3*(Q75-Q25),
  #         upper_thresh = Q75 + 3*(Q75-Q25)) %>% 
  #  select(Abbrev, lower_thresh, upper_thresh)
  
  #randPoints <- min(sum_stats$n, maxPoints)
  
  #rand_sample <- input_df %>%
  #  sample_n(randPoints) %>% 
  #  left_join(outlier_thresholds, by = 'Abbrev') %>% 
  #  mutate(outlier = ifelse(z<=upper_thresh & z>=lower_thresh, NA, "outlier")) %>% 
  #  select(Abbrev, z, outlier)
  
  if (i == 1) {
    all_sumstats <- sum_stats
    # all_randsamps <- rand_sample
  } else {
    all_sumstats <- bind_rows(all_sumstats, sum_stats)
    # all_randsamps <- bind_rows(all_randsamps, rand_sample)
  }
    
  setTxtProgressBar(pb, i)
  
}

write_csv(all_sumstats, "All-HUC8-Zstar-Summ-Stats.csv")

close(pb)

geo_files <- "F:/z-star-spatial-data/temp/z-star-uncertainty/"

# Get list of geotifs in the specified folder where the geospatial files live
file_list <- list.files(geo_files)
tif_list <- file_list[grepl(".tif", file_list) & !grepl(".tif.aux.xml", file_list) &  !grepl(".tif.xml", file_list)]

pb <- txtProgressBar(min = 0, max = length(tif_list), style = 3)

for (i in 1:length(tif_list)) {
  
  input_tif <- tif_list[i]
  
  print(input_tif)
  
  input_path <- paste(geo_files, input_tif, sep="")
  input_raster <- raster(input_path)
  input_df <- as.data.frame(rasterToPoints(input_raster)) 
  
  names(input_df) <- c("x", "y", "z")
  
  input_df <- input_df %>% mutate(Abbrev = str_remove(input_tif, ".tif"))
  
  sum_stats <- input_df %>%
    group_by(Abbrev) %>%
    summarise(n = n(), median = median(z))
  
  if (i == 1) {
    all_sumstats <- sum_stats
  } else {
    all_sumstats <- bind_rows(all_sumstats, sum_stats)
  }
  
  setTxtProgressBar(pb, i)
  
}

close(pb)

write_csv(all_sumstats, "All-HUC8-ZstarUncertainty-Medians.csv")


#remove_xtreme_outliers <- all_randsamps %>% filter(is.na(outlier))

#ggplot(data = all_randsamps, aes(x = Abbrev, y=z)) +
#  geom_boxplot(aes(color = as.factor(Abbrev))) +
#  geom_hline(data=data.frame(datum=c("MSL", "MHW"), z=c(0,1)), aes(yintercept = z, lty=datum))
  #geom_rug(aes(y=z, color = outlier))
