##----------------------
## Code to make combinations of models to test
## CM: 28/03/20
##
##----------------------

## likelihood code
like_code <- data.frame(dist = c("poisson", "negbin"),
                        code = c("y[i] ~ dpois(lambda[i])",
                                 "y[i] ~ dnegbin(p[i],r) \n p[i] <- r/(r+lambda[i])"
                                 ),
                        stringsAsFactors = FALSE)

ppi_code <- data.frame(dist = c("poisson", "negbin"),
                       code = c("res[i] <- y[i] - lambda[i] \n y.new[i] ~ dpois(lambda[i]) \n res.new[i] <- y.new[i] - lambda[i]",
                                "res[i] <- y[i] - lambda[i] \n  y.new[i] ~ dnegbin(p[i],r) \n res.new[i] <- y.new[i] - lambda[i]"
                                ),
                        stringsAsFactors = FALSE)


ppsum_code <- " fit <- sum(res[(T-5):T]) \n fit.new <- sum(res.new[(T-5):T])"

dist_prior_code <- data.frame(dist = c("poisson", "negbin"),
                              ##code = c("", "r ~ dunif(0,50)"),
                              code = c("", "r ~ dgamma(0.001, 0.001)"),
                              stringsAsFactors = FALSE)

## predictive distribution code (forecast)
pred_code <- data.frame(dist = c("poisson", "negbin"),
                        code = c("ypred[i] ~ dpois(lambda[T+i])",
                                 "p[T+i] <- r/(r+lambda[T+i]) \n ypred[i] ~ dnegbin(p[T+i], r)"
                                 ),
                        stringsAsFactors = FALSE)

## need to change the index here also
level_proc_code <- data.frame(level = c("const", "rw1", "rw2", "ar1"),
                              code = c(
                                  ## const
                                  "log_lambda[i] <- log_lambda[i-1] + delta[i-1]",
                                  ## rw1
                                  "log_lambda[i] ~ dnorm(log_lambda[i-1] + delta[i-1], tau_eps)",
                                  ## rw2
                                  "log_lambda[i] ~ dnorm(2 * log_lambda[i-1] - log_lambda[i-2] + delta[i-1], tau_eps)",
                                  ## ar1
                                  "log_lambda[i] ~ dnorm(phil * log_lambda[i-1] + delta[i-1], tau_eps)"
                              ),
                              stringsAsFactors = FALSE
                              )

## need to change the index here also
trend_proc_code <- data.frame(trend = c("const", "rw1", "rw2", "ar1"),
                              code = c(
                                  ## const
                                  "delta[i] <- delta[i-1]",
                                  ## rw1
                                  "delta[i] ~ dnorm(delta[i-1], tau_eta)",
                                  ## rw2
                                  "delta[i] ~ dnorm(2 * delta[i-1] - delta[i-2], tau_eta)",
                                  ## ar1
                                  "delta[i] ~ dnorm(phid * delta[i-1], tau_eta)"
                              ),
                              stringsAsFactors = FALSE
                              )

level_prior_code <- data.frame(level = c("const", "rw1", "rw2", "ar1"),
                               code = c(
                                   ## const
                                   "",
                                   c(rep(
                                       "sigma_eps ~ dunif(0, 1) \n tau_eps <- 1 / pow(sigma_eps, 2)", 2),
                                       "sigma_eps ~ dunif(0, 1) \n tau_eps <- 1 / pow(sigma_eps, 2) \n phil ~ dunif(-1, 1)")
                               ),
                              stringsAsFactors = FALSE)

trend_prior_code <- data.frame(trend = c("const", "rw1", "rw2", "ar1"),
                               code = c(
                                   ## const
                                   "",
                                   c(rep(
                                       "sigma_eta ~ dunif(0, 1) \n tau_eta <- 1 / pow(sigma_eta, 2)", 2),
                                       "sigma_eta ~ dunif(0, 1) \n tau_eta <- 1 / pow(sigma_eta, 2) \n phid ~ dunif(-1, 1)")
                               ),
                              stringsAsFactors = FALSE)

make_model_string <- function(dist, level, trend){
    ## likelihood string
    like_string <- paste(
        "for(i in 1:T) {",
        like_code[match(dist, like_code$dist), "code"],
        "lambda[i] <- exp(log_lambda[i])",
        ppi_code[match(dist, ppi_code$dist), "code"],
        "}",
        sep = "\n")
    ## predictive distribution string
    pred_string <- paste(
        "for(i in 1:Tf) {",
        pred_code[match(dist, pred_code$dist), "code"],
        "lambda[T+i] <- exp(log_lambda[T+i])",
        "}",
        sep = "\n")
    ## level process
    level_proc_string <- paste(
        ifelse(level == "rw2", "for(i in 3:(T+Tf)){", "for(i in 2:(T+Tf)){"),
        level_proc_code[match(level, level_proc_code$level), "code"],
        "}",
        sep = "\n")
    ## level process
    trend_proc_string <- paste(
        ifelse(trend == "rw2", "for(i in 3:(T+Tf)){", "for(i in 2:(T+Tf)){"),
        trend_proc_code[match(trend, trend_proc_code$trend), "code"],
        "}",
        sep = "\n")
    ## generic diffuse priors
    prior_generic_string <- "log_lambda[1] ~ dnorm(0, 10^-6) \ndelta[1] ~ dnorm(0, 10^-6)"
    if(level == "rw2"){
        prior_generic_string <- paste(prior_generic_string, "\n log_lambda[2] ~ dnorm(0, 10^-6)", sep = "\n")
    }
    if(trend == "rw2"){
        prior_generic_string <- paste(prior_generic_string, "\n delta[2] ~ dnorm(0, 10^-6)", sep = "\n")
    }    
    ## level prior code
    level_prior_string <- level_prior_code[match(level, level_prior_code$level), "code"]
    ## trend prior code
    trend_prior_string <- trend_prior_code[match(trend, trend_prior_code$trend), "code"]
    ## distribution prior code (negbin)
    dist_prior_string <- dist_prior_code[match(dist, dist_prior_code$dist), "code"]
    ##
    model_string <- paste(
        "model {\n",
        like_string,
        pred_string,
        level_proc_string,
        trend_proc_string,
        prior_generic_string,
        level_prior_string,
        trend_prior_string,
        dist_prior_string,
        ppsum_code,
        "}",
        sep = "\n")
    return(model_string)
}

