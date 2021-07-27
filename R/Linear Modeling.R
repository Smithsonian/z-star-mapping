# Linear modeling 

library(MuMIn)
library(tidyverse)
library(rgdal)
library(tmap)
library(sp)
library(gstat)
library(sjstats)
library(gridExtra)

# Spatial Trends
# Load spatial data we will use
z_star_sum_stats <- read_csv("tables/z-star-huc8-sum_stats.csv")

# Z star uncertainty 
z_star_uncertainty <- read_csv("data/z-star-watershed-summ-stats/All-HUC8-ZstarUncertainty-Medians.csv") %>% 
  rename(zStar_uncertainty = median) %>% 
  select(-n)

# Load up tidal range, and RSLR
mhw_msl <- read_csv("data/z-star-watershed-summ-stats/All-HUC8-MHW-MSL-Medians.csv") %>%
  rename(mhw_msl = median) %>% 
  select(-n)

rslr <- read_csv("data/z-star-watershed-summ-stats/All-HUC8-RSLR-Medians.csv") %>%
  rename(rslr = median) %>% 
  select(-n)

# Outliers were manually reviewed
outliers_reviewed <- read_csv("tables/Z-star-HUC8-outliers-manually-investigated.csv") %>% 
  filter(manual_check == "exclude")

full_modeling_file <- z_star_sum_stats %>% 
  left_join(z_star_uncertainty) %>% 
  left_join(mhw_msl) %>% 
  left_join(rslr) %>% 
  filter(complete.cases(.),
         ! Abbrev %in% outliers_reviewed$Abbrev) %>% 
  # Remove outliers
  mutate(IQR = Q75-Q25,
         log_MHW_MSL = log(mhw_msl),
         log_IQR= log(IQR),
         Coast = as.factor(Coast),
         log_ZstarUncertainty = log(zStar_uncertainty))

write_csv(full_modeling_file, "tables/Z-star-HUC8-full-modeling-file.csv")

# Visually check for normality
hist(full_modeling_file$log_ZstarUncertainty)
hist(full_modeling_file$log_IQR)
hist(full_modeling_file$median)
hist(full_modeling_file$rslr)
hist(full_modeling_file$log_MHW_MSL)

options(na.action = na.fail)
model1 <- lm(median ~ log_MHW_MSL*rslr, data = full_modeling_file)
model1_dredge <- dredge(model1)  

par(mar = c(3,5,6,4))
plot(model1_dredge, labAsExpr = TRUE)

# SpatialPointsDataFrame(points = data.frame(full_modeling_file$INSIDE_X,
#                                            full_modeling_file$INSIDE_Y))

best_model1 <- get.models(model1_dredge, 1)[[1]]
summary(best_model1)

best_model_Table1 <- anova_stats(best_model1)

best_model_Table1 <- best_model_Table1[order(-best_model_Table1$cohens.f),]

full_modeling_file$model1residuals <- best_model1$residuals
full_modeling_file$model1_fitted_values <- best_model1$fitted.values

ggplot(full_modeling_file, aes(x=INSIDE_X, y=INSIDE_Y)) +
  geom_point(aes(fill=model1residuals), pch=21, alpha=0.6) +
  scale_fill_gradientn(colours = rainbow(5))

ggplot(full_modeling_file, aes(x=exp(log_MHW_MSL), y=model1_fitted_values)) +
  geom_point(aes(color=rslr, pch=Coast), alpha=0.6) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(full_modeling_file, aes(x=exp(log_MHW_MSL), y=model1residuals)) +
  geom_point(aes(color=rslr, pch=Coast), alpha=0.6) +
  scale_color_gradientn(colours = rainbow(5))

model2 <- lm(log_IQR ~ log_MHW_MSL*rslr, data = full_modeling_file)
model2_dredge <- dredge(model2)  
head(model2_dredge)
plot(model2_dredge)
best_model2 <- get.models(model2_dredge, 1)[[1]]
best_model_Table2 <- summary(best_model2)

best_model2 <- get.models(model2_dredge, 1)[[1]]
summary(best_model2)

best_model_Table2 <- anova_stats(best_model2)

best_model_Table2 <- best_model_Table2[order(-best_model_Table2$cohens.f),]

full_modeling_file$model2residuals <- best_model2$residuals
full_modeling_file$model2_fitted_values <- best_model2$fitted.values

ggplot(full_modeling_file, aes(x=INSIDE_X, y=INSIDE_Y)) +
  geom_point(aes(fill=model2residuals, size=model2residuals), pch=21, alpha=0.6) +
  scale_fill_gradientn(colours = rainbow(5))

ggplot(full_modeling_file, aes(x=exp(log_MHW_MSL), y=exp(model2_fitted_values))) +
  geom_point(aes(color=rslr, pch=Coast), alpha=0.6) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(full_modeling_file, aes(x=exp(log_MHW_MSL), y=model2residuals)) +
  geom_point(aes(color=rslr, pch=Coast), alpha=0.6) +
  scale_color_gradientn(colours = rainbow(5))

model3 <- lm(log_ZstarUncertainty ~ log_MHW_MSL, data = full_modeling_file)
summary(model3)

best_model_Table3 <- anova_stats(model3)

best_model_Table3 <- best_model_Table3[order(-best_model_Table3$cohens.f),]

full_modeling_file$model3residuals <- model3$residuals
full_modeling_file$model3_fitted_values <- model3$fitted.values

ggplot(full_modeling_file, aes(x=log_MHW_MSL, y=log_ZstarUncertainty)) +
  geom_point()

###### Variograms for the three different model residuals ######## 
# Define the 1st order polynomial equation
coordinates(full_modeling_file) = ~INSIDE_X+INSIDE_Y
full_modeling_file@proj4string <- CRS("+proj=longlat +datum=WGS84")

full_modeling_file$X <- coordinates(full_modeling_file)[,1]
full_modeling_file$Y <- coordinates(full_modeling_file)[,2]

# Define the 1st order polynomial equation
f.1 <- as.formula(model1residuals ~ 1) 

# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
var.smpl.1 <- variogram(f.1, full_modeling_file)

# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit.1  <- fit.variogram(var.smpl.1, vgm("Sph"))

# The following plot allows us to assess the fit
plot(var.smpl.1, dat.fit.1, main = "Median Z* Residuals")

f.2 <- as.formula(model2residuals ~ 1) 

# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
var.smpl.2 <- variogram(f.2, full_modeling_file)

# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit.2  <- fit.variogram(var.smpl.2, vgm("Sph"))

# The following plot allows us to assess the fit
plot(var.smpl.2, dat.fit.2, main = "IQR Z* Residuals")

f.3 <- as.formula(model3residuals ~ 1) 

# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
var.smpl.3 <- variogram(f.3, full_modeling_file)

# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit.3  <- fit.variogram(var.smpl.3, vgm("Sph"))

# The following plot allows us to assess the fit
plot(var.smpl.3, dat.fit.3, main = "Z* Uncertainty Residuals")

####### Bootstrapping R2 values #####

# For each model, we are going to interate through the watersheds
tests_to_run <- c("model1residuals",
                  "model2residuals",
                  "model3residuals")
out_names <- c("Z* Median",
               "Z* Variability",
               "Z* Uncertainty")

store_outputs <- data.frame(test=rep(NA,length(tests_to_run)*nrow(full_modeling_file)),
                            HUC8=rep(NA,length(tests_to_run)*nrow(full_modeling_file)),
                            residual=rep(NA,length(tests_to_run)*nrow(full_modeling_file)),
                            krige.prediction=rep(NA,length(tests_to_run)*nrow(full_modeling_file)),
                            boot.mean=rep(NA,length(tests_to_run)*nrow(full_modeling_file)))
out_row = 1
for (j in 1:length(tests_to_run)) {
  
  for (i in 1:nrow(full_modeling_file)) {
    temp_leave_out <- full_modeling_file[i,]
    temp_include <- full_modeling_file[-i,]
    
    temp.f <- as.formula(eval(parse(text=tests_to_run[j])) ~ 1) 
    
    temp_semi_v <- variogram(temp.f, temp_include)
    temp_fit <- fit.variogram(temp_semi_v, vgm("Sph"))
    # plot(temp_semi_v, temp_fit)
    
    dat.krg <- krige(temp.f, temp_include, temp_leave_out, temp_fit)
    
    store_outputs$test[out_row] <- out_names[j]
    store_outputs$HUC8[out_row] <- full_modeling_file$Name[i]  
    store_outputs$residual[out_row] <- unlist(temp_leave_out@data[1, tests_to_run[j]])[1]
    store_outputs$krige.prediction[out_row] <- dat.krg$var1.pred
    store_outputs$boot.mean[out_row] <- mean(unlist(
      temp_include@data[, tests_to_run[j]]), na.rm=T)
    out_row = out_row + 1
  }
}

# We will leave out one point
# We will fit a semi-variogram using the rest of the points
# We will store the prediction
# At the end we will summarise the sum square errors, 
# Divide it by the variance (sum square xi - x-mean)
# Boot strapped "Pseudo-R2" is 1-(SSE/Var)

store_outputs_R2 <- store_outputs %>% 
  mutate(residual_minus_modeled_sq = (residual-krige.prediction)^2,
         residual_minus_bootmean_sq = (residual-boot.mean)^2) %>% 
  group_by(test) %>% 
  summarise(sse = sum(residual_minus_modeled_sq),
            var = sum(residual_minus_bootmean_sq)) %>% 
  mutate(R2 = 1 - (sse/var))

out_names <- c("Z* Median",
               "Z* Variability",
               "Z* Uncertainty")

best_model_Table1$test <- "Z* Median"
best_model_Table2$test <- "Z* Variability"
best_model_Table3$test <- "Z* Uncertainty"

best_covariate_models <- best_model_Table1 %>% 
  bind_rows(best_model_Table2) %>% 
  bind_rows(best_model_Table3) %>% 
  select(test, term, omegasq)

adjR2 <- best_covariate_models %>% 
  group_by(test) %>% 
  summarise(adjR2 = sum(omegasq, na.rm = T))

spatialR2 <- store_outputs_R2 %>% 
  left_join(adjR2) %>% 
  mutate(pseudo_R2 = (1-adjR2) * R2) %>% 
  select(test, pseudo_R2)

adjR2 <- adjR2 %>% mutate(term="Total covariate model",
                          stat = "adjusted-R^2",
                          model = "covariate") %>% 
  rename(variance_explained = adjR2)

spatialR2 <- spatialR2 %>% mutate(term="Total spatial model",
                    stat = "pseudo-R^2",
                    model = "spatial") %>% 
  rename(variance_explained = pseudo_R2)

best_covariate_models <- best_covariate_models %>% 
  mutate(term=str_replace(term, "log_MHW_MSL", "log(MHW-MSL)"),
         stat = "omega^2",
         model = "covariate") %>% 
  rename(variance_explained = omegasq) %>% 
  filter(complete.cases(.))

all_variance_explained <- best_covariate_models %>% 
  bind_rows(adjR2) %>% 
  bind_rows(spatialR2)

all_variance_explained$test <- factor(all_variance_explained$test, levels = out_names)

all_variance_explained$test2 <- factor(str_replace(all_variance_explained$test, "Z*\\*", "Z*'*' [MHW]* ' ' *"),
                                       levels = c("Z*'*' [MHW]* ' ' * Median", 
                                                  "Z*'*' [MHW]* ' ' * Variability", 
                                                  "Z*'*' [MHW]* ' ' * Uncertainty"))

ggplot(data = all_variance_explained, aes(x=term, y=variance_explained, color = model)) +
  geom_segment(aes(xend = term, yend = 0)) +
  geom_point() +
  facet_wrap(.~test2, labeller = "label_parsed") +
  ylab("Variance Explained") +
  xlab(element_blank()) +
  theme(legend.position = "none",
        axis.text.x=element_text(angle=45,hjust=1))

ggsave("figures/Modeling Effect sizes.pdf", width = 5, heigh = 3, unit = "in")
ggsave("figures/Modeling Effect sizes.jpg", width = 5, heigh = 3, unit = "in")


# median_real_data_plot
# IQR_real_data_plot
# zStarUncertaintyRealData

full_modeling_file_plots <- as.data.frame(full_modeling_file)

median_real_data_plot <- ggplot(full_modeling_file_plots, aes(x=mhw_msl, y=median)) +
  geom_point(aes(pch=Coast), size= 3, fill = 'darkgrey') +
  ylab(expression(paste("Z*"["MHW"], " Median", sep = ""))) +
  scale_shape_manual(values=c(21, 22, 23)) +
  xlab(NULL) +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))

(median_real_data_plot)

IQR_real_data_plot <- ggplot(full_modeling_file_plots, aes(x=mhw_msl, y=IQR)) +
  geom_point(aes(pch=Coast), size= 3, fill = 'darkgrey') +
  ylab(expression(paste("Z*"["MHW"], " IQR (log)", sep = ""))) +
  scale_shape_manual(values=c(21, 22, 23)) +
  scale_x_log10() +
  scale_y_log10() +
  xlab(NULL) +
  theme_minimal() +
  theme(legend.position = "none") 

(IQR_real_data_plot)

zStarUncertaintyRealData <- ggplot(full_modeling_file_plots, aes(x=mhw_msl, y=zStar_uncertainty)) +
  geom_point(aes(pch=Coast), size= 3, fill = 'darkgrey') +
  ylab(expression(paste("Z*"["MHW"], " Uncertainty (log)", sep = ""))) +
  scale_shape_manual(values=c(21, 22, 23)) +
  xlab("MHW-MSL (log)") +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(pch = element_blank())

(zStarUncertaintyRealData)

p <- arrangeGrob(median_real_data_plot, 
             IQR_real_data_plot,
             zStarUncertaintyRealData,
             ncol = 1,
             nrow = 3,
             heights = c(1.5,1,1))

ggsave("figures/Z* HUC8 Regressions.pdf", p, width = 3.54, height = 6, units = "in")
ggsave("figures/Z* HUC8 Regressions.jpg", p, width = 3.54, height = 6, units = "in")
