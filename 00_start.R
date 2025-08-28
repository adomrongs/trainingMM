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

# -------------------------------------------------------------------------
# Least-squares -----------------------------------------------------------
# -------------------------------------------------------------------------

n <- 12
n_blks <- 3
n_gens <- 4

# Design matrices
ff <- yield ~ 1 + block + gen
str(m <- model.frame(ff, data))
X <- model.matrix(ff, m)
y <- matrix(data[, "yield"])
print(X)
print(y)

# X = [1 X_blks X_gens]
X_intercept <- X[, 1, drop = FALSE]
X_blks <- X[, 2:3, drop = FALSE]
X_gens <- X[, 4:6, drop = FALSE]

# Betas
Xty <- t(X) %*% y
XtX <- t(X) %*% X
qr(XtX)$rank
XtX_inv <- solve(XtX)
beta <- XtX_inv %*% Xty
beta

# Number of effects
n_coef <- length(beta)
n_coef

# Sigma
y_hat <- X %*% beta
errors <- y - y_hat
SSE <- sum(errors^2)
SSE / (n - n_coef)

# R(beta) = beta'X'y
R_beta <- t(beta) %*% t(X) %*% y
R_mu <- t(beta[1, ]) %*% t(X_intercept) %*% y
R_blks <- t(beta[2:3, ]) %*% t(X_blks) %*% y
R_gens <- t(beta[4:6, ]) %*% t(X_gens) %*% y

SC_blks <- R_beta - R_mu - R_gens
SC_gens <- R_beta - R_mu - R_blks
SC_blks


# Using 'lm'
mod_0 <- lm(formula = yield ~ 1 + block + gen, data = data)
mod_0
anova(mod_0)

