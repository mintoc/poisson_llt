##---------------
## Take a look at covid-19 case growth rate
## Using count models with time-varying ln rate and trend
## predict from all models here
## CM: 25/3/2020
##---------------
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(R2jags)
##library(jagsUI)
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
##----------
## HOLD-BACK
##----------

nhold <- 5
sfit <- s[-((nrow(s)-(nhold-1)):nrow(s)),]
spred <- s[(nrow(s)-(nhold-1)):nrow(s),]

## response
n <- sfit$ncases
n[is.na(n)] <- 0 ## true zero or not?

T <- length(n)

## number of years to forecast (required in the general code)
Tf <- nhold

data <- list(y = n,
             T = T,
             Tf = Tf
             )

pars0 <- c("lambda", "delta", "sigma_eps", "sigma_eta", "ypred", "r", "phil", "phid", "fit", "fit.new")

all_fits_hold_df <- all_pred_hold_df <- pp_prob_df <- NULL

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
        ##plot(m)
        lambda_ci <- sapply(1:T, FUN = function(z){
            quantile(m$BUGSoutput$sims.array[, , paste0("lambda[", z, "]")], p = c(0.025, 0.975))
        })
        res <- sfit
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
        all_fits_hold_df <- rbind(all_fits_hold_df, res)
        ## forecasts
        ## prediction
        ypred_ci <- sapply(1:Tf, FUN = function(z){
            quantile(m$BUGSoutput$sims.array[, , paste0("ypred[", z, "]")], p = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95))
        })        
        ypred_median <- sapply(1:Tf, FUN = function(z){
            median(m$BUGSoutput$sims.array[, , paste0("ypred[", z, "]")])
        })
        ypred_df <- data.frame(date = sfit$date[T] + 1:Tf,
                               ypred_mean = m$BUGSoutput$mean$ypred,
                               ypred_median = ypred_median,
                               ypred_lwr50 = ypred_ci["25%",],
                               ypred_upr50 = ypred_ci["75%",],
                               ypred_lwr80 = ypred_ci["10%",],
                               ypred_upr80 = ypred_ci["90%",],
                               ypred_lwr90 = ypred_ci["5%",],
                               ypred_upr90 = ypred_ci["95%",],
                               model = mname)
        all_pred_hold_df <- rbind(all_pred_hold_df, ypred_df)
        ## hold back probabilities
        pp_prob <- sapply(1:nhold, function(z){
            samps <- m$BUGSoutput$sims.array[, , paste0("ypred[", z, "]")]
            obs <- spred$ncases[z]
            prob <- sum(samps >= 0.9 * obs & samps <= 1.1 * obs) / prod(dim(samps))
            return(prob)
        })
        hold_df <- data.frame(model = mname, date = spred$date, pp_prob = pp_prob)
        pp_prob_df <- rbind(pp_prob_df, hold_df)
        rm(list = c("res", "ypred_df", "hold_df"))
    }
}

## for side-by side vis and increasing complexity
all_fits_hold_df$fmodel <- factor(all_fits_hold_df$model,
                                  levels = paste(mods[,1], mods[,2], mods[,3], sep = "_"))

all_pred_hold_df$fmodel <- factor(all_pred_hold_df$model,
                             levels = paste(mods[,1], mods[,2], mods[,3], sep = "_"))

mnames <- paste(mods[,1], mods[,2], mods[,3], sep = "_")

pp_prob_df$fmodel <- factor(pp_prob_df$model,
                            levels = rev(mnames[order(mnames)]))

p3 <-
    ggplot(all_fits_hold_df, aes(x = date, y = ncases)) +
    geom_point(size = 2, pch = 1) +
    xlab("") +
    geom_ribbon(aes(ymin = lambda_lwr, ymax = lambda_upr), fill = "#F8766D", alpha = 0.5) +
    geom_line(aes(y = lambda)) +
    facet_wrap(~ fmodel, scales = "free") +
    ylab(paste0("Number of new cases in ", countryi)) +
    theme(strip.text.x = element_text(size = 12)) +
    geom_ribbon(data = all_pred_hold_df, aes(y = ypred_median, ymin = ypred_lwr90, ymax = ypred_upr90),
                fill = "thistle1") +
    geom_ribbon(data = all_pred_hold_df, aes(y = ypred_median, ymin = ypred_lwr80, ymax = ypred_upr80),
                fill = "plum2") +
    geom_ribbon(data = all_pred_hold_df, aes(y = ypred_median, ymin = ypred_lwr50, ymax = ypred_upr50),
                fill = "plum3") +
    geom_line(data = all_pred_hold_df, aes(y = ypred_mean, group = model), col = "plum4")+
    geom_line(data = all_pred_hold_df, aes(y = ypred_median, group = model), col = "plum4", linetype = 2) +
    geom_point(data = spred, pch = 8)

pdf(paste0("../tex/figures/all_fits_hold_", countryi,".pdf"), height = 9.6, width = 12.8)
print(p3)
dev.off()

pdf(paste0("../tex/figures/pp_prob_hold_", countryi,".pdf"), height = 8, width = 6)
ggplot(pp_prob_df, aes(x = date, y = fmodel, fill = pp_prob)) +
    geom_tile(, width=0.9, height=0.9) +
    scale_fill_gradient(low = "white", high = "red") +
    xlab("") + ylab("")
dev.off()

## difference from the mean/median
mdiff_df <- merge(all_pred_hold_df, spred[, c("ncases", "date")], by = "date")

## get the relative error 
mdiff_df$mean_re <- with(mdiff_df, 100 * (ypred_mean - ncases)/ncases)
mdiff_df$median_re <- with(mdiff_df, 100 * (ypred_median - ncases)/ncases)

library(reshape)
mdiff_long_df <- melt(mdiff_df[, c("date", "fmodel", "mean_re", "median_re")], id.vars = c("date", "fmodel"))

## cut off those more than 100% out
mdiff_long_df$value[abs(mdiff_long_df$value) > 100] <- NA
mdiff_long_df$label <- paste0(round(mdiff_long_df$value, 2), "%")
mdiff_long_df$label[is.na(mdiff_long_df$value)] <- "abs > 100%"

mdiff_long_df$fmodel <- factor(as.character(mdiff_long_df$fmodel),
                               levels = rev(mnames[order(mnames)]))

pdf(paste0("../tex/figures/mdiff_hold_", countryi,".pdf"), height = 10, width = 12)
ggplot(mdiff_long_df, aes(x = date, y = fmodel, fill = value)) +
    geom_tile(, width=0.9, height=0.9) +
    geom_text(aes(label = label)) +
    facet_wrap(~ variable) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_text(size=rel(1.7)),
          axis.text.y = element_text(size=rel(1.7)),
          strip.text.x = element_text(size = rel(1.7)))
dev.off()
