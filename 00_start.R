library(tidyverse)
library(emmeans)
library(lme4)

# RCBD
# 4 gens
# 3 blocks
# 12 observations
data <- read.csv("data/example_1.csv") |>
  mutate(gen = as.factor(gen), block = as.factor(block))
head(data)
str(data)

# Visualization
data |>
  ggplot(aes(x = gen, y = yield, color = block)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15)

# Mean by genotype
data |>
  group_by(gen) |>
  summarise(m = mean(yield))

# Mean by block
data |>
  group_by(block) |>
  summarise(m = mean(yield))

# Overall mean
data |>
  summarise(m = mean(yield))

# -------------------------------------------------------------------------
# Linear-Models -----------------------------------------------------------
# -------------------------------------------------------------------------

n <- 12
n_blks <- 3
n_gens <- 4

# Intercept ---------------------------------------------------------------

ff <- yield ~ 1
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas (mu)
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta_mu <- XtX_inv %*% Xty
beta_mu
y_hat <- X %*% beta_mu
errors <- y - y_hat
SSE_mu <- sum(errors^2)
SSE_mu

# Block -------------------------------------------------------------------

ff <- yield ~ -1 + block
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas (block)
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta_blk <- XtX_inv %*% Xty
beta_blk
y_hat <- X %*% beta_blk
errors <- y - y_hat
SSE_blk <- sum(errors^2)
SSE_blk

# Genotype ----------------------------------------------------------------

ff <- yield ~ -1 + gen
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas (genotype)
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta_gen <- XtX_inv %*% Xty
beta_gen
y_hat <- X %*% beta_gen
errors <- y - y_hat
SSE_gen <- sum(errors^2)
SSE_gen

# All effects -------------------------------------------------------------

ff <- yield ~ 1 + block + gen
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# Betas
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta <- XtX_inv %*% Xty
beta
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
SSE

# Reduction in SS if a factor is added
SSE_mu
SSE_mu - SSE_blk
SSE_mu - SSE_gen
SSE

# Sigma
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
sigma_2 <- SSE / (n - n_coef)
sigma_2

# Using 'lm'
mod <- lm(formula = yield ~ 1 + block + gen, data = data)
mod
anova(mod)
vcov(mod)
XtX_inv * sigma_2

# SST
SST <- t(y - beta_mu[1]) %*% (y - beta_mu[1])
SST
# -------------------------------------------------------------------------

# Number of effects
n_coef <- length(beta)
n_coef

# Overall mean
beta[1] + sum(c(0, beta[2:3])) / 3 + sum(c(0, beta[4:6])) / 4
beta_mu

# What about the gen missing level?
print(beta)
gens <- c("geng1" = 0, beta[4:6, ])
gens <- beta[1] + sum(c(0, beta[2:3])) / 3 + gens
gens

# What about the block missing level?
print(beta)
blks <- c("block1" = 0, beta[2:3, ])
blks <- beta[1] + sum(c(0, beta[4:6])) / 4 + blks
blks

# From the model
emmeans(mod, ~gen)
sqrt(XtX_inv[4:6, 4:6] * sigma_2)

emmeans(mod, pairwise ~ gen)

# -------------------------------------------------------------------------
# asreml ------------------------------------------------------------------
# -------------------------------------------------------------------------

asreml.options(Cfixed = TRUE)
mod_asr <- asreml(fixed = y ~ block + gen, data = data)
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
