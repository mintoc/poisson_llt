##---------------
## Take a look at covid-19 case growth rate
## Using count models with time-varying ln rate and trend
## fit all models here
## CM: 25/3/2020
##---------------
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(R2jags)
## library(jagsUI)
## data
d <- read.csv("https://raw.githubusercontent.com/mathsnuig/coronaviz/master/data/corona_island.csv")

## choose a country
##countryi <- "ireland"
countryi <- "italy"

s <- d %>%
    filter(country == countryi) %>%
    mutate(date = as.Date(date,format = "%d/%m/%Y")) %>%
    complete(date = seq.Date(min(date), max(date), by="day")) %>% ## get the missing days
    group_by(date) %>%
    summarise(ncases = sum(ncase), Date = min(date))

## response
n <- s$ncases
n[is.na(n)] <- 0 ## true zero or not?

T <- length(n)

## number of years to forecast (just used for plotting here)
Tf <- 3

## model
source("make_model_string.R")

## possible models
mods <- expand.grid(dist = c("poisson", "negbin"),
                    ##level = c("const", "rw1", "rw2", "ar1"),
                    level = c("const", "rw1", "rw2"),
                    ##trend = c("const", "rw1", "rw2", "ar1"),
                    trend = c("const", "rw1", "rw2"),
                    stringsAsFactors = FALSE
                    ##DIC = NA
)

data <- list(y = n,
             T = T,
             Tf = Tf
             )
pars0 <- c("lambda", "delta", "sigma_eps", "sigma_eta", "ypred", "r", "phil", "phid", "fit", "fit.new")

all_fits_df <- all_pred_df <- NULL

for(i in 1:nrow(mods)){
    print(i)
    ##i <- 5
    model_string <- make_model_string(dist = mods$dist[i],
                                      level = mods$level[i],
                                      trend = mods$trend[i])
    cat(model_string, file = "model.txt")
    m <- try(jags(data = data, n.iter=1e4, n.chains = 3, model.file = "model.txt",
                  parameters.to.save = pars0), silent = TRUE)
    if(class(m) != "try-error"){
        plot(m)
        lambda_ci <- sapply(1:T, FUN = function(z){
            quantile(m$BUGSoutput$sims.array[, , paste0("lambda[", z, "]")], p = c(0.025, 0.975))
        })
        res <- s
        res$lambda <- m$BUGSoutput$mean$lambda[1:T]
        res$lambda_lwr <- lambda_ci[1,]
        res$lambda_upr <- lambda_ci[2,]
        ## delta
        delta_mean <- m$BUGSoutput$mean$delta[1:T]
        delta_ci <- sapply(1:T, FUN = function(z){
            quantile(m$BUGSoutput$sims.array[, , paste0("delta[", z, "]")], p = c(0.025, 0.975))
        })
        res$delta <- delta_mean
        res$delta_lwr <- delta_ci[1,]
        res$delta_upr <- delta_ci[2,]
        ##
        res$DIC <- m$BUGSoutput$DIC
        mname <- paste(as.character(mods[i, ]), collapse = "_")
        res$model <- mname
        all_fits_df <- rbind(all_fits_df, res)
        ## forecasts
        ## prediction
        ypred_ci <- sapply(1:Tf, FUN = function(z){
            quantile(m$BUGSoutput$sims.array[, , paste0("ypred[", z, "]")], p = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95))
        })        
        ypred_median <- sapply(1:Tf, FUN = function(z){
            median(m$BUGSoutput$sims.array[, , paste0("ypred[", z, "]")])
        })
        ypred_df <- data.frame(date = s$date[T] + 1:Tf,
                               ypred_mean = m$BUGSoutput$mean$ypred,
                               ypred_median = ypred_median,
                               ypred_lwr50 = ypred_ci["25%",],
                               ypred_upr50 = ypred_ci["75%",],
                               ypred_lwr80 = ypred_ci["10%",],
                               ypred_upr80 = ypred_ci["90%",],
                               ypred_lwr90 = ypred_ci["5%",],
                               ypred_upr90 = ypred_ci["95%",],
                               model = mname)
        all_pred_df <- rbind(all_pred_df, ypred_df)
        rm(list = c("res", "ypred_df"))
    }
}

## for side-by side vis and increasing complexity
all_fits_df$fmodel <- factor(all_fits_df$model,
                             levels = paste(mods[,1], mods[,2], mods[,3], sep = "_"))

dic_df <- subset(all_fits_df, date == date[1])
dic_df$ncases <- Inf

dic_df$dDIC <- with(dic_df, paste("dDIC:", round(DIC - min(DIC), 2)))

p0 <-
    ggplot(all_fits_df, aes(x = date, y = ncases)) +
    geom_point(size = 2, pch = 1) +
    xlab("") +
    geom_ribbon(aes(ymin = lambda_lwr, ymax = lambda_upr), fill = "#F8766D", alpha = 0.5) +
    geom_line(aes(y = lambda)) +
    facet_wrap(~ fmodel, scales = "free") +
    geom_text(data = dic_df, aes(label = dDIC), hjust = 0, vjust = 2) +
    ylab(paste0("Number of new cases in ", countryi)) +
    theme(strip.text.x = element_text(size = 12))
    ##ggtitle("Posterior Poisson rate (lambda) and 95% credible intervals")

pdf(paste0("../tex/figures/all_fits_", countryi,".pdf"), height = 9.6, width = 12.8)
print(p0)
dev.off()

## plot of the deltas
p1 <- ggplot(all_fits_df, aes(x = date, y = delta)) +
    geom_line(colour = "darkblue") +
    xlab("") +
    geom_ribbon(aes(ymin = delta_lwr, ymax = delta_upr), fill = "#00BFC4", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_wrap(~ fmodel) +
    geom_text(data = dic_df, aes(label = dDIC, y = 0.6), hjust = 0) +
    ylab("Posterior log rate change (delta)") +
    theme(strip.text.x = element_text(size = 12))
    

pdf(paste0("../tex/figures/all_deltas_", countryi,".pdf"), height = 9.6, width = 12.8)
p1
dev.off()

## plot of the forecast
all_pred_df$fmodel <- factor(all_pred_df$model,
                             levels = paste(mods[,1], mods[,2], mods[,3], sep = "_"))

p2 <-
    p0 +
    geom_ribbon(data = all_pred_df, aes(y = ypred_median, ymin = ypred_lwr90, ymax = ypred_upr90),
                fill = "thistle1") +
    geom_ribbon(data = all_pred_df, aes(y = ypred_median, ymin = ypred_lwr80, ymax = ypred_upr80),
                fill = "plum2") +
    geom_ribbon(data = all_pred_df, aes(y = ypred_median, ymin = ypred_lwr50, ymax = ypred_upr50),
                fill = "plum3") +
    geom_line(data = all_pred_df, aes(y = ypred_mean, group = model), col = "plum4")+
    geom_line(data = all_pred_df, aes(y = ypred_median, group = model), col = "plum4", linetype = 2)

pdf(paste0("../tex/figures/all_pred_", countryi,".pdf"), height = 9.6, width = 12.8)
p2
dev.off()
