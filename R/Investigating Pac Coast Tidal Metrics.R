# Testing the hypothesis that HAT*, HOT*, and DHQ* scale with tidal amplitude

library(tidyverse)

othet_pac_stats <- read_csv("data/HAT DHQ Gauges Used in Janousek et al 2019 Reanalysis.csv")

othet_pac_stats_star <- othet_pac_stats %>% 
  mutate(HAT_star = (HAT-MSL)/(MHW-MSL),
         HOT_star = (HOT-MSL)/(MHW-MSL),
         DHQ_star = (DHQ-MSL)/(MHW-MSL),
         log_tidal_amp = log(MHW-MSL),
         log_HAT_star = log(HAT_star),
         log_HOT_star = log(HOT_star),
         log_DHQ_star = log(DHQ_star))

model1 <- lm(log_HAT_star~log_tidal_amp, data = othet_pac_stats_star)
summary(model1)

model2 <- lm(log_HOT_star~log_tidal_amp, data = othet_pac_stats_star)
summary(model2)

model3 <- lm(log_DHQ_star~log_tidal_amp, data = othet_pac_stats_star)
summary(model3)

plot_1 <- ggplot(othet_pac_stats_star, aes(x = MHW, y = HAT_star)) +
  geom_point(pch = 21, size = 4, fill = "white") +
  xlab("MHW-MSL (log scale)") +
  scale_x_log10() +
  geom_smooth(method = "lm", se = T, color = "white", lty = 2) +
  ylab(expression(paste("HAT*"["MHW"], " (log scale)"))) +
  theme_dark() +
  ggtitle("A.")
  # geom_label(data = data.frame(x = 1, y = 2.2), aes(x=x, y=y), label = expression(paste("p = 0.0003, R"^2, " = 0.81", sep="")))

plot_1  

plot_2 <- ggplot(othet_pac_stats_star, aes(x = MHW, y = HOT_star)) +
  geom_point(pch = 21, size = 4, fill = "white") +
  xlab("MHW-MSL (log scale)") +
  scale_x_log10() +
  geom_smooth(method = "lm", se = T, color = "white", lty = 2) +
  ylab(expression(paste("HOT*"["MHW"], " (log scale)"))) +
  theme_dark() +
  ggtitle("B.")
  # geom_label(data = data.frame(x = 1, y = 2.66), aes(x=x, y=y), label = expression(paste("p = 0.034, R"^2, " = 0.38", sep="")))

plot_2


plot_3 <- ggplot(othet_pac_stats_star, aes(x = MHW, y = DHQ_star)) +
  geom_point(pch = 21, size = 4, fill = "white") +
  xlab("MHW-MSL (log scale)") +
  scale_x_log10() +
  geom_smooth(method = "lm", se = T, color = "white", lty = 2) +
  ylab(expression(paste("DHQ*"["MHW"], " (log scale)"))) +
  theme_dark() +
  ggtitle("C.")
  # geom_label(data = data.frame(x = 1, y = 0.3), aes(x=x, y=y), label = expression(paste("p < 0.0001, R"^2, " = 0.86", sep="")))

plot_3

gridExtra::grid.arrange(plot_1, plot_2, plot_3, ncol = 3)

g <- gridExtra::arrangeGrob(plot_1, plot_2, plot_3, ncol = 3)

ggsave(file="figures/Supplemental_Fig1.pdf", g,
       width = 7.25,
       height = 2.5) #saves g

ggsave(file="figures/Supplemental_Fig1.jpg", g,
       width = 7.25,
       height = 2.5) #saves g
