
library(foreign)
library(tidyverse)
library(gridExtra)

{
  # Function that returns Unbiased Root Mean Squared Error
  # From Jolliff et al., 2009
  standardizedSigma <-function(referenceValues, modeledValues) {
    sigma_r <- sd(referenceValues)
    sigma_m <- sd(modeledValues)
    return(sigma_m / sigma_r)
  }
  
  correlationCoef <- function(referenceValues, modeledValues) {
    n <- length(modeledValues)
    mu_m <- mean(modeledValues)
    mu_r <- mean(referenceValues)
    sigma_m <- sd(modeledValues)
    sigma_r <- sd(referenceValues)
    R <- ((sum((modeledValues - mu_m) * (referenceValues - mu_r)) 
           / n) 
          / (sigma_m * sigma_r))
    return(R)
  }
  
  unbiasedRmse <- function(referenceValues,modeledValues) { 
    r<-referenceValues
    m<-modeledValues
    m_bar <- mean(modeledValues) # turned off NA'rms so you have to make the two vectors equal ahead of time
    r_bar <- mean(referenceValues)
    rmse <- sqrt(mean( ((m-m_bar) - (r - r_bar) )^2))
    sigma_m <- sd(modeledValues)
    sigma_r <- sd(referenceValues)
    rmse <- rmse * sign(sigma_m - sigma_r)
    return(rmse)
  }
  
  standardizedUnbiasedRmse <- function(referenceValues, modeledValues) {
    sigma_ast <- standardizedSigma(referenceValues, modeledValues)
    R <- correlationCoef(referenceValues, modeledValues)
    sigma_m <- sd(modeledValues)
    sigma_r <- sd(referenceValues)
    unbiasedRMSE_ast <- sqrt(1 + sigma_ast^2 - (2 * sigma_ast * R))
    unbiasedRMSE_ast <- unbiasedRMSE_ast * sign(sigma_m - sigma_r)
    return(unbiasedRMSE_ast)
  }
  
  bias <- function(referenceValues, modeledValues) { return(mean(modeledValues)-mean(referenceValues)) }
  
  standardizedBias <- function(referenceValues, modeledValues) {
    sigma_r <- sd(referenceValues)
    return(bias(referenceValues, modeledValues) / sigma_r)
  }
  
}

original <- as.tibble(
  read.dbf("data/independentValidation/GCREW-TMON-indy-validation-Original-LiDAR.dbf", as.is=T)) %>%
  rename(referenceZ = Elevation, modeledZ.original = RASTERVALU)

lean <- as.tibble(
  read.dbf("data/independentValidation/GCREW-TMON-indy-validation-Lean-Corrected.dbf", as.is=T)) %>%
  rename(referenceZ = Elevation, modeledZ.lean = RASTERVALU)

vegClasses <- as.tibble(
  read.dbf("data/independentValidation/GCREW-TMON-indy-validation-VegClass-Corrected.dbf", as.is=T)) %>%
  rename(referenceZ = Elevation, modeledZ.vegClasses = RASTERVALU)

krigged <- as.tibble(
  read.dbf("data/independentValidation/GCREW-TMON-indy-validation-KriggedDem.dbf", as.is=T)) %>%
  rename(referenceZ = Elevation, modeledZ.krigged = RASTERVALU)

compiledStrategies <- original %>%
  left_join(lean) %>%
  left_join(vegClasses) %>%
  left_join(krigged) %>%
  mutate(modeledZ.simple = modeledZ.original - 0.173) %>%
  select(-Code)

compiledStrategies[compiledStrategies == -9999] <- NA
compiledStrategies[compiledStrategies == "NaN"] <- NA 

compiledStrategies <- compiledStrategies %>%  
  filter(complete.cases(.))

plot(compiledStrategies$referenceZ, compiledStrategies$modeledZ.original)
plot(compiledStrategies$referenceZ, compiledStrategies$modeledZ.lean)
plot(compiledStrategies$referenceZ, compiledStrategies$modeledZ.vegClasses)
plot(compiledStrategies$referenceZ, compiledStrategies$modeledZ.krigged)
plot(compiledStrategies$referenceZ, compiledStrategies$modeledZ.simple)

outputPlot = data.frame(analysis = c("Simple Offset", "LEAN", "Veg. Class", "Krig"),
                        bias = rep(NA,4),
                        unbiasedRmse = rep(NA, 4),
                        standardizedBias = rep(NA, 4),
                        standardizedUnbiasedRmse = rep(NA, 4))


# Simple Strategy

outputPlot$bias[1] <- bias(referenceValues=compiledStrategies$referenceZ, 
                           modeledValues=compiledStrategies$modeledZ.simple)

outputPlot$unbiasedRmse[1] <- unbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                           modeledValues=compiledStrategies$modeledZ.simple)

outputPlot$standardizedBias[1] <- standardizedBias(referenceValues=compiledStrategies$referenceZ, 
                                                   modeledValues=compiledStrategies$modeledZ.simple)

outputPlot$standardizedUnbiasedRmse[1] <- standardizedUnbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                                   modeledValues=compiledStrategies$modeledZ.simple)

# Lean Strategy


outputPlot$bias[2] <- bias(referenceValues=compiledStrategies$referenceZ, 
                           modeledValues=compiledStrategies$modeledZ.lean)

outputPlot$unbiasedRmse[2] <- unbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                           modeledValues=compiledStrategies$modeledZ.lean)

outputPlot$standardizedBias[2] <- standardizedBias(referenceValues=compiledStrategies$referenceZ, 
                                                   modeledValues=compiledStrategies$modeledZ.lean)

outputPlot$standardizedUnbiasedRmse[2] <- standardizedUnbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                                                   modeledValues=compiledStrategies$modeledZ.lean)

# Veg Class Strategy 

outputPlot$bias[3] <- bias(referenceValues=compiledStrategies$referenceZ, 
                           modeledValues=compiledStrategies$modeledZ.vegClasses)

outputPlot$unbiasedRmse[3] <- unbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                           modeledValues=compiledStrategies$modeledZ.vegClasses)

outputPlot$standardizedBias[3] <- standardizedBias(referenceValues=compiledStrategies$referenceZ, 
                                                   modeledValues=compiledStrategies$modeledZ.vegClasses)

outputPlot$standardizedUnbiasedRmse[3] <- standardizedUnbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                                                   modeledValues=compiledStrategies$modeledZ.vegClasses)

# Veg Class Strategy 

outputPlot$bias[4] <- bias(referenceValues=compiledStrategies$referenceZ, 
                           modeledValues=compiledStrategies$modeledZ.krigged)

outputPlot$unbiasedRmse[4] <- unbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                           modeledValues=compiledStrategies$modeledZ.krigged)

outputPlot$standardizedBias[4] <- standardizedBias(referenceValues=compiledStrategies$referenceZ, 
                                                   modeledValues=compiledStrategies$modeledZ.krigged)

outputPlot$standardizedUnbiasedRmse[4] <- standardizedUnbiasedRmse(referenceValues=compiledStrategies$referenceZ, 
                                                                   modeledValues=compiledStrategies$modeledZ.krigged)


total_max = max(abs(c(outputPlot$standardizedUnbiasedRmse, outputPlot$standardizedBias)))


ciclesWeWant <- c(0.25, 0.5, 0.75, 1)
circleLWD <- c("<1", "<1", "<1", "1+")
for (i in 1:length(ciclesWeWant)) {
  circle1pts.x1 <- seq(0, ciclesWeWant[i], by=0.01)
  circle1ptx.y1 <- sqrt(ciclesWeWant[i]^2 - circle1pts.x1^2)
  tempCircleDF = data.frame(x=c(circle1pts.x1, rev(circle1pts.x1), -circle1pts.x1, -rev(circle1pts.x1)), y=c(circle1ptx.y1, -rev(circle1ptx.y1), -circle1ptx.y1, rev(circle1ptx.y1)))
  tempCircleDF["circleDefinition"] <- rep(toString(ciclesWeWant[i]), nrow(tempCircleDF))
  tempCircleDF["circleLWD"] <- rep(toString(circleLWD[i]), nrow(tempCircleDF))
  if (i == 1) {
    circleDF = tempCircleDF
  } else {
    circleDF<-rbind(circleDF, tempCircleDF)
  }
}

outputPlot <- outputPlot %>%
  mutate(standardizedTotalRmse = sqrt(standardizedBias^2 + standardizedUnbiasedRmse^2))

#palette using grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#palette using black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


targetDiagram <- ggplot(data = outputPlot, aes(x = standardizedUnbiasedRmse, y = standardizedBias)) +
  geom_polygon(data=circleDF, color="grey", fill=NA, aes(x=x, y=y, mapping=circleDefinition, lwd=circleLWD)) +
  geom_point(aes(color = analysis, shape = analysis), size=4) +
  guides(col=guide_legend(title=NULL), shape = guide_legend(title=NULL)) +
  scale_x_continuous(limits=c(-total_max, total_max)) +
  scale_y_continuous(limits=c(-total_max, total_max)) + 
  xlab(expression(paste("Root Mean Square Error'* "%*%" sign(", sigma["m"], "-", sigma["r"], ")"))) +
  ylab("Bias*") +
  scale_colour_manual(values = cbbPalette) +
  scale_size_manual(values=c(.25,1), guide=FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("A. Precision vs Accuracy")

outputPlot["Scalability"] <- 1:4

tradeoffs <- ggplot(outputPlot, aes(x=Scalability, y=standardizedTotalRmse)) +
  geom_point(stat="identity", size=4, aes(color = analysis, pch = analysis)) + 
  geom_segment(aes(y=standardizedTotalRmse, 
                   yend=0, 
                   x=Scalability, 
                   xend=Scalability, color=analysis)) + 
  xlab("Most Scalable to Least") +
  ylab("Performance (Root Mean Square Error*)") +
  theme_bw() +
  guides(col=guide_legend(title=NULL), shape = guide_legend(title=NULL)) +
  scale_colour_manual(values = cbbPalette) +
  theme(axis.text.x=element_blank()) +
  ggtitle("B. Scalability vs Performance")
  
grid.arrange(targetDiagram, tradeoffs, nrow=1)
