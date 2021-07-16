# Plots for Z* and Z* uncertainty with HUC8's

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(sf)
library(rgdal)
library(gridExtra)
library(tmaptools)
library(MuMIn)

## Set up stuff for whole map #######

# Load spatial data we will use
huc8shp <- st_read("data/z-star-watershed-summ-stats/HUC8_Tidal_170404.shp")
huc3aea <- st_transform(huc8shp, crs = "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0")
b <- st_bbox(huc3aea)

# Load tabular data we will use
z_star <- read_csv("data/z-star-watershed-summ-stats/All-HUC8-Zstar-Summ-Stats.csv")

z_star_uncertainty <- read_csv("data/z-star-watershed-summ-stats/All-HUC8-ZstarUncertainty-Medians.csv")

# join huc8 data to spatial data
huc8_aea <- huc3aea %>% 
  select(Name, HUC8, Abbrev, States, Coast, INSIDE_X, INSIDE_Y, geometry)

# Filter into three coast maps.
# One coast at a time 

# Arrange West coast from north to south
pacific <- huc8_aea %>% 
  filter(Coast == "Pacific") %>% 
  arrange(-INSIDE_Y) %>% 
  # Create integers
  mutate(int_position = 1:length(INSIDE_X))

gc_start <- max(pacific$int_position) + 1

# Arrange Gulf coast from east to west
gulf <- huc8_aea %>% 
  filter(Coast == "Gulf") %>%
  arrange(INSIDE_X) %>%
  # Create integers
  mutate(int_position = gc_start:(gc_start+length(INSIDE_X)-1))

ac_start <- max(gulf$int_position) + 1

# Arrange atlantic coast from south to north
atlantic <-  huc8_aea %>% 
  filter(Coast == "Atlantic") %>% 
  arrange(INSIDE_Y) %>% 
  # Create integers
  mutate(int_position = ac_start:(ac_start+length(INSIDE_X)-1))

# Stitch files back together
coasts_in_color_order <- rbind(pacific, gulf, atlantic) %>%
  # Order them in integer order
  arrange(int_position)

huc8_info <- data.frame(Abbrev = coasts_in_color_order$Abbrev,
                        Name = coasts_in_color_order$Name,
                        States = coasts_in_color_order$States,
                        Coast = coasts_in_color_order$Coast,
                        INSIDE_X = coasts_in_color_order$INSIDE_X,
                        INSIDE_Y = coasts_in_color_order$INSIDE_Y,
                        int_position = coasts_in_color_order$int_position)

# Look at the color brewer pallates
display.brewer.all()

# Use color brewer to assign spectral palletes
coasts_in_color_order["man_colors"] <- get_brewer_pal(palette="Spectral", 
                                                      n = max(coasts_in_color_order$int_position),
                                                      stretch = TRUE)

map.sf <- ne_countries(scale = 'medium', type = 'map_units',
                       returnclass = 'sf')

map.na.sf <- map.sf[map.sf$continent == 'North America',]

map.na.sf.aea <- st_transform(map.na.sf, 
                              crs = "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0")

center_map_figure <- ggplot(data = coasts_in_color_order) + 
  geom_sf(data=map.na.sf.aea, color="black", size=0.1, fill="white") +
  geom_sf(fill=coasts_in_color_order$man_colors, color="black", size=0.05) +
  coord_sf(xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]),
           crs="+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0") +
  # Probably want to take off x and y axes text
  theme_dark() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# ggsave("center_map.pdf")




## Let's try to remove outliers ######

hist(z_star$median)

z_star_sum_stats <- z_star %>% 
  group_by() %>%
  summarise(median_median = median(median), 
            median_Q25 = quantile(median, 0.25),
            median_Q75 = quantile(median, 0.75),
            median_upper_outlier = median_Q75 + 3*(median_Q75-median_Q25),
            median_lower_outler = median_Q25 - 3*(median_Q75-median_Q25),
            median_min = min(median),
            median_max = max(median),
            median_Q975 = median(Q975),
            median_Q75 = median(Q75),
            median_Q25 = median(Q25),
            median_Q025 = median(Q025)
  ) 

z_star_sum_stats_view <- huc8_info %>% 
  left_join(z_star) %>%
  arrange(-n) %>% 
  mutate(pct_of_all_data = n / sum(z_star$n) * 100)

View(z_star_sum_stats_view %>% arrange(int_position))  

write_csv(z_star_sum_stats_view %>% arrange(int_position),
          "tables/z-star-huc8-sum_stats.csv")

nrow(z_star %>% filter(median >= 1)) / nrow(z_star)

z_star_outlier_analysis <- z_star %>% 
  group_by() %>%
  mutate(median_Q25 = quantile(median, 0.25),
         median_Q75 = quantile(median, 0.75),
         upper_outlier = median_Q75 + 3*(median_Q75-median_Q25),
         lower_outler = median_Q25 - 3*(median_Q75-median_Q25),
         outlier = ifelse(median >= lower_outler & median <= upper_outlier, "in range", "outlier"))

big_differences <- z_star_outlier_analysis %>%
  group_by(outlier) %>%
  summarise(min_z_star = min(median),
            max_z_star = max(median),
            mean_z_star = mean(median),
            min_n = min(n),
            max_n = max(n),
            mean_n = mean(n))

z_star_outliers <- z_star_outlier_analysis %>% filter(outlier == "outlier") %>% 
  select(Abbrev)

z_star_outlier_view <- z_star_outlier_analysis %>% 
  filter(outlier == "outlier") %>% 
  left_join(huc8_info)

write_csv(z_star_outlier_view, "tables/Z-star-HUC8-outliers.csv")

# Outliers were manually reviewed
outliers_reviewed <- read_csv("tables/Z-star-HUC8-outliers-manually-investigated.csv") %>% 
  filter(manual_check == "exclude")

z_star <- z_star %>% filter(!Abbrev %in% outliers_reviewed$Abbrev)
z_star_uncertainty <- z_star_uncertainty %>% filter(!Abbrev %in% outliers_reviewed$Abbrev)

## Z star uncertainty ############

# Join uncertainty table with shapefile
graph_zstar_sigma <- coasts_in_color_order %>% 
  left_join(z_star_uncertainty) %>% 
  filter(!Abbrev %in% outliers_reviewed$Abbrev)

# May want to consider removing outlier watersheds
# Do these watersheds have somehting in common

# Create a lollipop graph for all plots for median uncertainty
ggplot(data = graph_zstar_sigma, aes(x=int_position, y=median)) +
  geom_point(pch=16, color = graph_zstar_sigma$man_colors) +
  geom_segment(aes(xend=int_position, yend=0), color = graph_zstar_sigma$man_colors) +
  ylab(expression(paste("Uncertainty (Z*"["MHW"],")",sep=""))) +
  geom_hline(data=data.frame(datum=c("MHW"), z_star=c(1)), aes(yintercept=z_star), lty=1) +
  # scale_y_log10() +
  theme_dark() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Each coast needs it's own variation on this plot
# With a filtered datasets data set
# east and west coast need to coordinates flipped
# We can try reverse the axes to make the plot look cool

# West

wc_uncertainty <- graph_zstar_sigma %>% 
  filter(Coast == "Pacific")

# Create a lollipop graph for all plots for median uncertainty
wc_uncertainty_plot <- ggplot(data = wc_uncertainty, aes(x=int_position, y=median)) +
  geom_point(pch=16, color = wc_uncertainty$man_colors) +
  geom_segment(aes(xend=int_position, yend=0), color = wc_uncertainty$man_colors) +
  ylab(NULL) +
  geom_hline(data=data.frame(datum=c("MHW"), z_star=c(1)), aes(yintercept=z_star), lty=1) +
  # scale_y_log10() +
  theme_dark() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() +
  scale_y_reverse() +
  scale_x_reverse()

(wc_uncertainty_plot)

# Gulf 
gc_uncertainty <- graph_zstar_sigma %>% 
  filter(Coast == "Gulf")

gc_uncertainty_plot <- ggplot(data = gc_uncertainty, aes(x=int_position, y=median)) +
  geom_point(pch=16, color = gc_uncertainty$man_colors) +
  geom_segment(aes(xend=int_position, yend=0), color = gc_uncertainty$man_colors) +
  ylab(expression(paste("Uncertainty (Z*"["MHW"],")",sep=""))) +
  geom_hline(data=data.frame(datum=c("MHW"), z_star=c(1)), aes(yintercept=z_star), lty=1) +
  # scale_y_log10() +
  theme_dark() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # coord_flip() +
  scale_y_reverse()

(gc_uncertainty_plot)

# East

ec_uncertainty <- graph_zstar_sigma %>% 
  filter(Coast == "Atlantic")

ec_uncertainty_plot <- ggplot(data = ec_uncertainty, aes(x=int_position, y=median)) +
  geom_point(pch=16, color = ec_uncertainty$man_colors) +
  geom_segment(aes(xend=int_position, yend=0), color = ec_uncertainty$man_colors) +
  ylab(NULL) +
  geom_hline(data=data.frame(datum=c("MHW"), z_star=c(1)), aes(yintercept=z_star), lty=1) +
  # scale_y_log10() +
  theme_dark() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip() 

(ec_uncertainty_plot)

# then we set up a grid plot
# If we envision 12 quadrats four columns and three rows

# 1 2 2 3
# 1 2 2 3
# NA 4 4 4

# 1 is pacific coast
# 2 is center map
# 3 is atlantic coast
# 4 is gulf coast
# 5 and 6 are NA plots

plot_m<-matrix(c(1,NA,2,4,3,NA),nrow=2, ncol=3)

allGrobs <- list(wc_uncertainty_plot, center_map_figure, ec_uncertainty_plot,
                 gc_uncertainty_plot)

allMaps <- arrangeGrob(grobs = allGrobs, layout_matrix = plot_m,
                       heights=c(1.5,1),
                       widths=c(1,1.5,1),
                       padding=unit(0.25, "line"))

ggsave("Z-star-uncertainty-HUC8-3-coast.pdf", allMaps, width = 7.25, height=3.5)                  
ggsave("Z-star-uncertainty-HUC8-3-coast.jpg", allMaps, width = 7.25, height=3.5)                  


## Z star Distributions ######

# Join uncertainty table with shapefile
graph_zstar <- coasts_in_color_order %>% 
  left_join(z_star, by="Abbrev") %>%
  # Graph the outer wiskers
  # It will either be Q1 or Q3 - 1.5 x IQR or min/max depending on which is less
  mutate(upper_outlier = Q75 + 1.5*(Q75-Q25),
         lower_outlier = Q25 - 1.5*(Q75-Q25),
         upper_wisker = NA,
         lower_wisker = NA,
         IQR = (Q75-Q25))

for(i in 1:nrow(graph_zstar)) {
  graph_zstar$upper_wisker[i] <- min(c(graph_zstar$upper_outlier[i], 
                                            graph_zstar$max[i]))
  graph_zstar$lower_wisker[i] <- max(c(graph_zstar$lower_outlier[i], 
                                             graph_zstar$min[i]))
}

# May want to consider removing outlier watersheds
# Do these watersheds have somehting in common

# Create a boxplot graph for all plotss
ggplot(data = graph_zstar, aes(x=int_position, 
                               lower = Q25, 
                               upper = Q75,
                               middle = median,
                               ymin = lower_wisker,
                               ymax = upper_wisker)) +
  geom_boxplot(color = graph_zstar$man_colors, stat = "identity") +
  ylab(expression(paste("Z*"["MHW"]))) +
  geom_hline(data=data.frame(datum=c("MHW", "MSL"), z_star=c(1, 0)), aes(yintercept=z_star, lty=datum)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(data = graph_zstar, aes(x=INSIDE_Y, 
                               y=median)) +
  geom_point(color = graph_zstar$man_colors, aes(shape=Coast)) +
  facet_wrap(.~Coast) +
  xlab("Latitude (dd)") +
  ylab(expression(paste("median (Z*"["MHW"], ")", sep="")))

ggplot(data = graph_zstar, aes(x=INSIDE_Y, 
                               y=IQR)) +
  geom_point(color = graph_zstar$man_colors, aes(shape=Coast)) +
  facet_wrap(.~Coast) +
  xlab("Latitude (dd)") +
  ylab(expression(paste("IQR (Z*"["MHW"], ")", sep="")))

# Each coast needs it's own variation on this plot
# With a filtered datasets data set
# east and west coast need to coordinates flipped
# We can try reverse the axes to make the plot look cool

# West

wc_zstar <- graph_zstar %>% 
  filter(Coast == "Pacific")

# Create a boxplot graph for all plots for median uncertainty
wc_zstar_plot <- ggplot(data = wc_zstar, aes(x=int_position, 
                               lower = Q25, 
                               upper = Q75,
                               middle = median,
                               ymin = lower_wisker,
                               ymax = upper_wisker)) +
  geom_boxplot(color = rev(wc_zstar$man_colors), stat = "identity", fill=NA,
               mapping = aes(group = int_position)) +
  ylab(expression(paste("Z*"["MHW"]," (Higher | Lower)",sep=""))) +
  geom_hline(data=data.frame(datum=c("MHW", "MSL"), z_star=c(1, 0)), aes(yintercept=z_star, lty=datum)) +
  theme_dark() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        # axis.title.x=element_blank(),
        legend.position = "left") +
  coord_flip() +
  scale_y_reverse() +
  scale_x_reverse() 
   
(wc_zstar_plot)

wc_bb <- st_bbox(wc_zstar)

wc_huc8 <- ggplot(data = wc_zstar) + 
  geom_sf(data=map.na.sf.aea, color="black", size=0.1, fill="white") +
  geom_sf(fill=wc_zstar$man_colors, color="black", size=0.05) +
  coord_sf(xlim = c(wc_bb["xmin"], wc_bb["xmax"]), ylim = c(wc_bb["ymin"], wc_bb["ymax"]),
           crs="+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0") +
  # Probably want to take off x and y axes text
  xlab(" ") +
  theme_dark() +
  theme(#axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.title.y=element_blank()
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank()
        )

wc_both_plots <- arrangeGrob(grobs = list(wc_zstar_plot, wc_huc8), nrow = 1, ncol=2,
            widths=c(2,1),
            padding=unit(0.25, "line"))

ggsave("WC_HUC8_plots.pdf", wc_both_plots, height=5, width=5)
ggsave("WC_HUC8_plots.jpg", wc_both_plots, height=5, width=5)


# Gulf 
gc_zstar <- graph_zstar %>% 
  filter(Coast == "Gulf")

gc_zstar_plot <- ggplot(data = gc_zstar, aes(x=int_position, 
                                             lower = Q25, 
                                             upper = Q75,
                                             middle = median,
                                             ymin = lower_wisker,
                                             ymax = upper_wisker)) +
  geom_boxplot(color = gc_zstar$man_colors, stat = "identity", fill = NA, mapping = aes(group = int_position)) +
  ylab(expression(paste("Z*"["MHW"]," (Higher | Lower)",sep=""))) +
  geom_hline(data=data.frame(datum=c("MHW", "MSL"), z_star=c(1, 0)), aes(yintercept=z_star, lty=datum)) +
  theme_dark() + 
  theme(# axis.title.y=element_blank(),
    axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")



(gc_zstar_plot)

gc_bb <- st_bbox(gc_zstar)

gc_huc8 <- ggplot(data = gc_zstar) + 
  geom_sf(data=map.na.sf.aea, color="black", size=0.1, fill="white") +
  geom_sf(fill=gc_zstar$man_colors, color="black", size=0.05) +
  coord_sf(xlim = c(gc_bb["xmin"], gc_bb["xmax"]), ylim = c(gc_bb["ymin"], gc_bb["ymax"]),
           crs="+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0") +
  # Probably want to take off x and y axes text
  xlab(" ") +
  theme_dark() +
  theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    axis.title.y=element_blank()
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank()
  )

gc_both_plots <- arrangeGrob(grobs = list(gc_huc8, gc_zstar_plot), nrow = 2, ncol=1,
                             heights = c(1,1.25),
                             padding=unit(0.25, "line"))

ggsave("GC_HUC8_plots.pdf", gc_both_plots, height=5, width=5)
ggsave("GC_HUC8_plots.jpg", gc_both_plots, height=5, width=5)


# East

ec_zstar <- graph_zstar %>% 
  filter(Coast == "Atlantic")

ec_zstar_plot <- ggplot(data = ec_zstar, aes(x=int_position, 
                                             lower = Q25, 
                                             upper = Q75,
                                             middle = median,
                                             ymin = lower_wisker,
                                             ymax = upper_wisker)) +
  geom_boxplot(color = ec_zstar$man_colors, stat = "identity", fill = NA,
               mapping=aes(group = int_position)) +
  ylab(expression(paste("Z*"["MHW"]," (Higher | Lower)",sep=""))) +
  geom_hline(data=data.frame(datum=c("MHW", "MSL"), z_star=c(1, 0)), aes(yintercept=z_star, lty=datum)) +
  theme_dark() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        # axis.title.x=element_blank(),
        legend.position = "left")  +
  coord_flip() + 
  scale_y_reverse()
  # scale_x_reverse() 

(ec_zstar_plot)

ec_bb <- st_bbox(ec_zstar)

head(arrange(ec_zstar, -upper_outlier))

ec_huc8 <- ggplot(data = ec_zstar) + 
  geom_sf(data=map.na.sf.aea, color="black", size=0.1, fill="white") +
  geom_sf(fill=ec_zstar$man_colors, color="black", size=0.05) +
  coord_sf(xlim = c(ec_bb["xmin"], ec_bb["xmax"]), ylim = c(ec_bb["ymin"], ec_bb["ymax"]),
           crs="+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0") +
  # Probably want to take off x and y axes text
  xlab(" ") +
  theme_dark() +
  theme(#axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    axis.title.y=element_blank()
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank()
  )

ec_both_plots <- arrangeGrob(grobs = list(ec_zstar_plot, ec_huc8), nrow = 1, ncol=2,
                             widths = c(1,1),
                             padding=unit(0.25, "line"))

ggsave("EC_HUC8_plots.pdf", ec_both_plots, height=5.5, width=5)
ggsave("EC_HUC8_plots.jpg", ec_both_plots, height=5.5, width=5)


# then we set up a grid plot
# If we envision 12 quadrats four columns and three rows

# 1 2 2 3
# 1 2 2 3
# 5 4 4 6

# 1 is pacific coast
# 2 is center map
# 3 is atlantic coast
# 4 is gulf coast
# 5 and 6 are NA plots

plot_m<-matrix(c(1,NA,2,4,3,NA),nrow=2, ncol=3)

allGrobs <- list(wc_zstar_plot, center_map_figure, ec_zstar_plot,
                 gc_zstar_plot)

allMaps <- arrangeGrob(grobs = allGrobs, layout_matrix = plot_m,
                       heights=c(2,1),
                       widths=c(1,2,1),
                       padding=unit(0.25, "line"))

ggsave("Z-star-HUC8-3-coast.pdf", allMaps, width = 7.25, height=3.5)                  

