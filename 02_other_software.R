library(tidyverse)
library(emmeans)
library(lme4)
library(MASS)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

# -------------------------------------------------------------------------
# asreml ------------------------------------------------------------------
# -------------------------------------------------------------------------

library(asreml)
asreml.options(Cfixed = TRUE)
mod_asr <- asreml(fixed = yield ~ block + gen, data = data)
summary(mod_asr)

# Coefficients
mod_asr$coefficients$fixed
mod_asr$Cfixed
preds <- predict(mod_asr, classify = "gen", vcov = TRUE, sed = TRUE)
preds$pvals
preds$vcov
preds$sed
sqrt(sigma_2 / 3) # SE of the BLUE
sqrt(sigma_2 / 3 + sigma_2 / 3) # SE of the difference
