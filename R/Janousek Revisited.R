# Janousek revisited 
library(tidyverse)
library(maps)
library(RColorBrewer)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# devtools::install_github("https://github.com/yutannihilation/ggsflabel")
library(ggsflabel)

all_data <- read_csv("data/Janoueseck_et_al_2019_plots.csv")

site_data <- all_data %>% 
  mutate(zStar = (elevation-MSL)/(MHW-MSL),
         site_id = recode(site_id, 
                          "Upper Newport Bay"="Newport Bay",
                          "Sweetwater Marsh"= "San Diego Bay")) %>% 
  select(site_id, plot_id, longitude, latitude, MHW, MSL, zStar) %>% 
  unique() %>% 
  group_by(site_id, plot_id) %>%
  summarise(longitude = first(longitude),
            latitude = first(latitude),
            MHW = first(MHW),
            MSL = first(MSL),
            zStar = first(zStar)) %>% 
  ungroup()

site_data_summary <- site_data %>% 
  group_by(site_id, MHW, MSL) %>% 
  summarise(median_Zstar = median(zStar),
            IQR = IQR(zStar),
            n = n(),
            longitude = median(longitude),
            latitude = median(latitude)) %>% 
  mutate(tidal_amp = MHW - MSL,
         log_tidal_amp = log(tidal_amp),
         log_IQR = log(IQR)) %>% 
  arrange(latitude)

model_1 <- lm(log_IQR ~ log_tidal_amp, data = site_data_summary)
summary(model_1)

model_2 <- lm(median_Zstar ~ log_tidal_amp, data = site_data_summary)
summary(model_2)

zStarDatums <- data.frame(zStar=c(0, 1),
                          datum=c("MSL", "MHW"))

all_group_latitudes <- site_data %>% 
  group_by(site_id) %>% 
  summarise(latitude = median(latitude)) %>% 
  arrange(-latitude)

all_data_plot <- site_data %>% 
  mutate(site_id = factor(site_id, all_group_latitudes$site_id))

site_data_summary_plot <- site_data_summary %>% 
  ungroup() %>% 
  mutate(site_id = factor(site_id, all_group_latitudes$site_id))

distFig <- ggplot(data=all_data_plot, aes(y=zStar, x=site_id)) +
  geom_boxplot(fill=NA, aes(color = site_id), show.legend = FALSE) +
  scale_color_brewer(palette = "Paired") +
  xlab("") + 
  # scale_y_continuous(limits = c(-0.5,5)) +
  geom_hline(data=zStarDatums, aes(yintercept = zStar, lty=datum)) +
  ylab(expression(paste("Z*"["MHW"]))) +
  xlab(NULL) +
  theme_dark() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom") +
  ggtitle("B.")
  
cor_median <- ggplot(site_data_summary_plot, aes(x = tidal_amp, y = median_Zstar)) +
  geom_point(aes(fill = site_id), show.legend = F, pch = 21, size = 4) +
  scale_fill_brewer(palette = "Paired") +
  xlab("MHW-MSL (log scale)") +
  scale_x_log10() +
  geom_smooth(method = "lm", se = T, color = "white", lty = 2) +
  ylab(expression(paste("Median Z*"["MHW"]))) +
  theme_dark() +
  ggtitle("C.") +
  geom_label(data = data.frame(x = 0.75, y = 2.2), aes(x=x, y=y), label = expression(paste("p = 0.0011, R"^2, " = 0.67", sep="")))
  
cor_median  
  
cor_IQR <- ggplot(site_data_summary_plot, aes(x = tidal_amp, y = IQR)) +
  geom_point(aes(fill = site_id), show.legend = F, pch = 21, size = 4) +
  scale_fill_brewer(palette = "Paired") +
  xlab("MHW-MSL (log scale)") +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = T, color = "white", lty = 2) + 
  ylab(expression(paste("IQR Z*"["MHW"], " (log scale)", sep=""))) +
  theme_dark() +
  ggtitle("D.") +
  geom_label(data = data.frame(x = 0.75, y = 0.825), aes(x=x, y=y), label = expression(paste("p = 0.044, R"^2, " = 0.35", sep="")))


cor_IQR  
  
(distFig)

library(maps)

sites <- all_data %>% 
  group_by(site_id) %>% 
  summarise(lat = median(latitude),
            long = median(longitude))

map.sf <- ne_countries(scale = 'medium', type = 'map_units',
                       returnclass = 'sf')

map.na.sf <- map.sf[map.sf$continent == 'North America',]

map.na.sf.aea <- st_transform(map.na.sf, 
                              crs = "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0")

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

site_data_summary_plot_sf <- st_as_sf(site_data_summary_plot,
                                           coords = c("longitude", "latitude"),
                                           crs = projcrs)

site_data_summary_plot_aea <- st_transform(site_data_summary_plot_sf,
                                                crs = "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0")

site_data_summary_plot_aea_labs <- site_data_summary_plot_aea %>% 
  filter(site_id %in% c("Tijuana Estuary", "San Diego Bay"))

b <- st_bbox(site_data_summary_plot_aea)

mapFig <- ggplot(data = site_data_summary_plot_aea) +
  geom_sf(data=map.na.sf.aea, color="black", size=0.1, fill="white") +
  geom_sf(aes(fill=site_id), color="black",
           show.legend = T, pch = 21, size = 4) +
  geom_sf_label_repel(aes(label = site_id, fill = site_id),
                      color = 'black',
                      # box.padding = NA,
                      point.padding = NA,
                      show.legend = FALSE,
                      size = 2,
                      force = 100) +
  scale_fill_brewer(palette = "Paired") +
  coord_sf(xlim = c(b["xmin"]-200000, b["xmax"]+80000), ylim = c(b["ymin"]-0.5, b["ymax"]+0.5),
           crs="+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0") +
  xlab(NULL) +
  ylab(NULL) +
  theme_dark() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  guides(color=guide_legend(ncol=2)) +
  ggtitle("A.") 

(mapFig)

gridExtra::grid.arrange(mapFig, distFig, cor_median, cor_IQR, nrow=2, ncol=2,
                        heights = c(3,2))

g <- gridExtra::arrangeGrob(mapFig, distFig, cor_median, cor_IQR, nrow=2, ncol=2,
                            heights = c(3,2))
ggsave(file="figures/Janousek 2016 Reanalysis.pdf", g,
       width = 7.25,
       height = 7.25) #saves g

ggsave(file="figures/Janousek 2016 Reanalysis.jpg", g,
       width = 7.25,
       height = 7.25) #saves g
