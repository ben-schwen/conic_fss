# install.packages("data.table")
# install.packages("MASS")
library(data.table)
require(MASS)

simulate_data <- function(n, beta, family = gaussian(), emu = 0, esd = 0, ...) {
    sigma = diag(length(beta) - 1L)
    family_name = family$family
    k <- length(beta) - 1L
    x <- MASS::mvrnorm(n, mu = double(k), Sigma = sigma)
    e <- rnorm(n, mean=emu, sd=esd)
    eta <- as.vector(cbind(1, x) %*% beta) + e
    mu <- family$linkinv(eta)
    simulate <- list(
        gaussian = function(n, mu, sd = 1, ...) rnorm(n = n, mean = mu, sd = sd),
        binomial = function(n, mu, ...) rbinom(n, 1, mu),
        poisson = function(n, mu, ...) rpois(n = n, lambda = mu)
    )[[family_name]]
    y <- simulate(n, mu)
    colnames(x) <- sprintf("x%s", seq_len(NCOL(x)))
    dat <- cbind(y = y, as.data.frame(x))
    attr(dat, "mu") <- mu
    attr(dat, "true_beta") <- beta
    dat
}

# gaussian
set.seed(1)
sim <- CJ(N=1000, SD=c(0.05, 0.10, 0.50, 1.00, 5.00), NVAR=seq(20, 50, 5), family="gaussian", RUN=1:5)
sim[, ID := .I]

data_gauss <- vector(mode="list", length=nrow(sim)+1)
data_gauss[[1]] <- sim
for (i in seq(nrow(sim))) {
    n <- sim[i,N]
    sd <- sim[i, SD]
    nvar <- sim[i, NVAR]
    k <- ceiling(nvar/2)
    uk <- nvar-k
    data <- simulate_data(n, c(1, rep(1,k), rep(0,uk)), esd=sd, family=gaussian())
    data <- data[c(seq(data)[-1], 1)]
    data_gauss[[i+1]] <- data
}
names(data_gauss) <- c("meta", sprintf("gauss%03d", seq(nrow(sim))))
dir.create('data', showWarnings = FALSE)
saveRDS(data_gauss, file="data/data_gauss.rds")

# binomial
set.seed(1)
sim <- CJ(N=1000, SD=c(0), NVAR=seq(5, 40, 5), family="binomial", RUN=1:5)
sim[, ID := .I]

data_bin <- vector(mode="list", length=nrow(sim)+1)
data_bin[[1]] <- sim
for (i in seq(nrow(sim))) {
    n <- sim[i,N]
    sd <- sim[i, SD]
    nvar <- sim[i, NVAR]
    k <- ceiling(nvar/2)
    uk <- nvar-k
    data <- simulate_data(n, c(1, rep(1,k), rep(0,uk)), esd=sd, family=binomial())
    data <- data[c(seq(data)[-1], 1)]
    data_bin[[i+1]] <- data
}
names(data_bin) <- c("meta", sprintf("binomial%03d", seq(nrow(sim))))
dir.create('data', showWarnings = FALSE)
saveRDS(data_bin, file="data/data_bin.rds")

# poisson
set.seed(1)
sim <- CJ(N=1000, SD=c(0), NVAR=seq(5, 40, 5), family="poisson", RUN=1:5)
sim[, ID := .I]

data_poisson <- vector(mode="list", length=nrow(sim)+1)
data_poisson[[1]] <- sim
for (i in seq(nrow(sim))) {
    n <- sim[i,N]
    sd <- sim[i, SD]
    nvar <- sim[i, NVAR]
    k <- ceiling(nvar/2)
    uk <- nvar-k
    data <- simulate_data(n, c(1, rep(1,k), rep(0,uk)), esd=sd, family=poisson())
    data <- data[c(seq(data)[-1], 1)]
    data_poisson[[i+1]] <- data
}
names(data_poisson) <- c("meta", sprintf("poisson%02d", seq(nrow(sim))))
dir.create('data', showWarnings = FALSE)
saveRDS(data_poisson, file="data/data_poisson.rds")

