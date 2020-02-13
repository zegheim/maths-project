###########################
# CONFIGURATION VARIABLES #
###########################

fname <- "~/Documents/diss/data/df.csv"
seed <- 17071996L
working.dir <- "~/Documents/diss/"

### DO NOT EDIT BELOW THIS LINE ###

###########
# IMPORTS #
###########

library(MASS)
library(optimx)
library(pscl)

########
# MAIN #
########

setwd(working.dir)
set.seed(seed)

Y <- df$count
# Normalise all covariates to have mean 0 and sd 1
X_scaled <- scale(as.matrix(subset(df, select = -c(count))))

# Add intercept column
X <- cbind(rep(1, length(Y)), X_scaled)
colnames(X)[1] <- "intercept"

loglik <- function(par, X, Y) {
  XB_PLUS <- X %*% par[1:5]
  XB_ZERO <- X %*% par[6:10]
  sum(
    ifelse(Y == 0, log(1 - inv.logit(XB_ZERO)), 0) +
      ifelse(Y != 0, log(inv.logit(XB_ZERO)) - exp(XB_PLUS) + (Y*XB_PLUS) - lgamma(Y + 1) - log(1 - exp(-exp(XB_PLUS))), 0)
  )                                    
}

par_init <- as.numeric(c(
  glm.fit(X, Y, family = poisson())$coefficients,
  glm.fit(X, factor(Y>0), family = binomial(link = "logit"))$coefficients
))

opt <- optimx(par_init, loglik, Y = Y, X = X, method = "BFGS", itnmax = 10000, control=list(maximize = TRUE))
B_PLUS_hat <- as.numeric(tail(opt, 1)[1:10])[1:5]
B_ZERO_hat <- as.numeric(tail(opt, 1)[1:10])[6:10]

# Check with using library
fire.hurdle <- hurdle(count ~ scale(x) + scale(y) + scale(elevation) + scale(max.temp), data = df, separate = FALSE)

# Compare final value
opt$value; fire.hurdle$optim$value

# Compare coefficients
as.numeric(fire.hurdle$coefficients$count); B_PLUS_hat
as.numeric(fire.hurdle$coefficients$zero); B_ZERO_hat

inf_matrix <- -1 * optimHess(par = as.numeric(opt[1:10]), fn = loglik, Y = Y, X = X, control=list(fnscale = -1))
inf_matrix[1:5, 6:10] <- 0
inf_matrix[6:10, 1:5] <- 0
inv.inf_matrix <- solve(inf_matrix)

# generating samples of MLE
mles <- mvrnorm(n = 100, mu = as.numeric(tail(opt, 1)[1:10]), Sigma = inv.inf_matrix)