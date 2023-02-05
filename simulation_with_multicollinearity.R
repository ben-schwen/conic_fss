.rsim_mcga <- function(n, p, alpha, standardize = TRUE) {
    mc_style <- TRUE
    if (mc_style) {
        p1 <- p + 1L
        Z <- matrix(rnorm(n * p1), n, p1)
    } else {
        p1 <- p
        Z <- matrix(rnorm(n * p1), n, p1)
    }
    X <- sqrt((1 - alpha^2)) * Z[, seq_len(p)] + alpha * Z[, p1]
    colnames(X) <- sprintf("A%i", seq_len(NCOL(X)))
    if (isTRUE(standardize)) {
        X <- scale(X)
    }
    ev <- eigen(t(X) %*% X)
    # beta_0 is taken to be zero
    beta <- ev$vectors[, which.max(ev$values)]
    list(X = X, beta = beta)
}


# This function reproduces the simulation setting from McDonald et al. (1975).
#
# @param n an integer giving the number of observations to be generated.
# @param alpha a double giving the parameter $\alpha$ which controls the pairwise correlations.
# @param sd a double giving the standard deviation to be used when generating Gaussian data.
#        This parameter is passed to the \code{rnorm} function and controls the standard deviation
#        of the error term $\epsilon$
# @param p an integer giving the number of covariates to be generated.
# @param standardize a logical, if \code{TRUE} the design matrix $X$ is standardized otherwise it is not.
#
# @references
# A Monte Carlo Evaluation of Some Ridge-Type Estimators
# Gary C. McDonald, Diane I. Galarneau
# Journal of the American Statistical Association
# Vol. 70, No. 350 (Jun., 1975), pp. 407-416 (10 pages)
# https://www.jstor.org/stable/2285832
# https://doi.org/10.2307/2285832
rsim_mcga <- function(n, alpha, sd = 0.1, family = gaussian(), p = 3L, standardize = TRUE) {
    dat <- .rsim_mcga(n = n, p = p, alpha = alpha, standardize = standardize)
    X <- dat[["X"]]
    beta <- dat[["beta"]]
    eta <- drop(X %*% beta)
    mu <- family$linkinv(eta)
    if (isTRUE(family$family == "gaussian")) {
        y <- rnorm(nrow(X), mu, sd = sd)
    } else if (isTRUE(family$family == "binomial")) {
        y <- rbinom(nrow(X), size = 1, prob = mu)
    } else if (isTRUE(family$family == "poisson")) {
        y <- rpois(nrow(X), lambda = mu)
    } else {
        stop("unkown family")
    }

    list(data = cbind(as.data.frame(X), y = y), beta = c(0, beta))
}



data_dir <- normalizePath("data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

nsim <- 1000L
nobs <- 100L
alphas <- c(0.7, 0.8, 0.9, 0.95, 0.99)

#
# Simulate gaussian data
#
family <- gaussian()
sigmas <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1)

params <- expand.grid(nobs = nobs, alpha = alphas, sd = sigmas, seed = seq_len(nsim),
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
dim(params)

for (i in seq_len(nrow(params))) {
    ofi <- file.path(data_dir, sprintf("%s_%05i.rds", family$family, i))
    if (file.exists(ofi)) {
        next()
    }
    param <- as.list(params[i, ])
    set.seed(param[["seed"]])
    data <- rsim_mcga(n = param[["nobs"]], alpha = param[["alpha"]], sd = param[["sd"]], family = family)
    saveRDS(data, ofi)
}


#
# Simulate binomial data
#
family <- binomial()

params <- expand.grid(nobs = nobs, alpha = alphas, seed = seq_len(nsim),
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
dim(params)

for (i in seq_len(nrow(params))) {
    ofi <- file.path(data_dir, sprintf("%s_%05i.rds", family$family, i))
    if (file.exists(ofi)) {
        next()
    }
    param <- as.list(params[i, ])
    set.seed(param[["seed"]])
    data <- rsim_mcga(n = param[["nobs"]], alpha = param[["alpha"]], family = family)
    saveRDS(data, ofi)
}


#
# Simulate poisson data
#
family <- poisson()

params <- expand.grid(nobs = nobs, alpha = alphas, seed = seq_len(nsim),
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
dim(params)

for (i in seq_len(nrow(params))) {
    ofi <- file.path(data_dir, sprintf("%s_%05i.rds", family$family, i))
    if (file.exists(ofi)) {
        next()
    }
    param <- as.list(params[i, ])
    set.seed(param[["seed"]])
    data <- rsim_mcga(n = param[["nobs"]], alpha = param[["alpha"]], family = family)
    saveRDS(data, ofi)
}
