library(Colossus)
library(data.table)
library(survival)
library(dplyr)
data(cancer, package = "survival")
veteran %>% setDT()
df <- copy(veteran)
# Make the same adjustments as Epicure example 6.5
karno <- df$karno
karno[93] <- 20
df$karno <- karno / 100
df$trt <- df$trt - 1
df$trt <- as.integer(df$trt == 0)
cell_lvl <- c("large", "squamous", "smallcell", "adeno")
df$cell <- as.integer(factor(df$celltype, level = cell_lvl)) - 1
df$karno50 <- df$karno - 50
control <- list(ncores = 1, maxiter = 20, halfmax = 1)
#
a_n <- c(-1, 1.17, -0.01)
model <- Pois_Strata(time, status, cell) ~ loglinear(trt, 0) + loglin - dose(karno, 1) + PA()
poisres <- PoisRun(model, df, a_n = a_n, control = control, keep_constant = c(0,1,0))
print(poisres$LogLik)
