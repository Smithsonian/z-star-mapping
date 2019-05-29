
#devtools::install_github("sgsokol/Deriv")

library(Deriv)

zstar <- function(E,MHW,MSL) { (MHW-E)/(MHW-MSL) }

dfdE <- Deriv(zstar, "E")

dfdMHW <- Deriv(zstar, "MHW", cache.exp = F)

dfdMSL <- Deriv(zstar, "MSL")

