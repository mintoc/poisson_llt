##---------------
## Take a look at covid-19 new case growth rate
## Using a Poisson with time-varying ln rate and trend
## CM: 25/3/2020
## Inspired by https://github.com/mintoc/environmetrics/blob/master/tmbIntro/dynamic_poisson.R and using modified JAGS code originally written by Andrew Jackson (TCD)
##---------------
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(R2jags)

## data
d <- read.csv("https://raw.githubusercontent.com/mathsnuig/coronaviz/master/data/corona_island.csv")

## sum the cases and add the zeros
s <- d %>%
    filter(country=="ireland") %>%
    mutate(date = as.Date(date,format = "%d/%m/%Y")) %>%
    complete(date = seq.Date(min(date), max(date), by="day")) %>% ## get the missing days
    group_by(date) %>%
    summarise(ncases = sum(ncase), Date = min(date))

## response
n <- s$ncases
n[is.na(n)] <- 0 ## true zero or not?

T <- length(n)

## JAGS model
modelstring <- "
model {
for(i in 1:N) { 
  y[i] ~ dpois(lambda[i])
  lambda[i] <- exp(log_lambda[i])
}
for(i in 2:N){
  log_lambda[i] ~ dnorm(log_lambda[i-1] + delta[i-1], tau_eps)
  delta[i] ~ dnorm(delta[i-1], tau_eta)
}
## diffuse prior on the first time point log_lambda[1] and delta[1]
log_lambda[1] ~ dnorm(0, 10^-6)
delta[1] ~ dnorm(0, 10^-6)
## prior on process error variances
sigma_eps ~ dunif(0, 1)
sigma_eta ~ dunif(0, 1)
##sigma_eps ~ dt(0, .04, 1)I(0,) ## sensitivity to variance priors
##sigma_eta ~ dt(0, .04, 1)I(0,)
tau_eps <- 1 / pow(sigma_eps, 2)
tau_eta <- 1 / pow(sigma_eta, 2)
}
"

data <- list(y = n,
             N = T
             )

cat(modelstring, file = "model.txt")

##
m <- jags(data = data, n.iter=1e5, n.chains = 3, model.file = "model.txt",
          parameters.to.save = c("lambda", "delta", "sigma_eps", "sigma_eta"))

plot(m)

## posteriors credible intervals
## lambda
lambda_ci <- sapply(1:T, FUN = function(z){
    quantile(m$BUGSoutput$sims.array[, , paste0("lambda[", z, "]")], p = c(0.025, 0.975))
})

s$lambda <- m$BUGSoutput$mean$lambda
s$lambda_lwr <- lambda_ci[1,]
s$lambda_upr <- lambda_ci[2,]

## delta
delta_mean <- m$BUGSoutput$mean$delta
delta_ci <- sapply(1:T, FUN = function(z){
    quantile(m$BUGSoutput$sims.array[, , paste0("delta[", z, "]")], p = c(0.025, 0.975))
})

s$delta <- delta_mean
s$delta_lwr <- delta_ci[1,]
s$delta_upr <- delta_ci[2,]

p0 <- ggplot(s, aes(x = date, y = ncases)) +
    geom_point(size = 2, pch = 1) +
    xlab("") +
    geom_ribbon(aes(ymin = lambda_lwr, ymax = lambda_upr), fill = "#F8766D", alpha = 0.5) +
    geom_line(aes(y = lambda)) +
    ylab("Number of new cases in Ireland") +
    ggtitle("Posterior Poisson rate (lambda) and 95% credible intervals")

p1 <- ggplot(s, aes(x = date, y = delta)) +
    geom_point(size = 2, colour = "cadetblue") +
    geom_linerange(aes(ymin = delta_lwr, ymax = delta_upr), colour = "cadetblue") +
    ggtitle("Posterior log rate change (delta)")

## most recent
deltaT_samp <- as.numeric(m$BUGSoutput$sims.array[, , paste0("delta[", T, "]")])
prop_increase <- exp(deltaT_samp) - 1
df <- data.frame(prop_increase)
prop_increase_mean <- mean(prop_increase)

p2 <- ggplot(df, aes(x = prop_increase)) +
    geom_density(fill = "skyblue", colour = NA) +
    geom_vline(xintercept = prop_increase_mean) +
    xlab("New case growth rate") +
    ggtitle(paste("New case growth rate", format(s$date[T], "%B %d %Y"), ":", round(prop_increase_mean, 3))) 

pdf("new_case_poisson_analysis.pdf", height = 8, width = 10)
plot(m)
p0; p1; p2
dev.off()

## probability that delta is decreasing (linear change in delta over time)
## slope over days for each sample
dims <- dim(m$BUGSoutput$sims.array)
nsamp <- dims[1]
nchains <- dims[2]

slopes <- matrix(NA, nrow = nsamp, ncol = nchains)

for(i in 1:nsamp){
    for(j in 1:nchains){
        y <- m$BUGSoutput$sims.array[i, j, grep("delta", dimnames(m$BUGSoutput$sims.array)[[3]])]
        X <- cbind(1, 1:T)
        beta <- solve(t(X) %*% X) %*% t(X) %*% matrix(y)
        slopes[i,j] <- beta[2,1]
    }
}

matplot(slopes, type = "l", lty = 1)

prop.table(table(slopes < 0))
