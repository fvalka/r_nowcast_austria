# remotes::install_github("joachim-gassen/tidycovid19")
library(tidyverse)
library(tidycovid19)
library(rstanarm)
library(ggplot2)
library(ggpubr)
library(xts)
library(gtools)
library(data.table)
library(lubridate)
library(EpiNow)
library(EpiSoon)
library(tictoc)
library(furrr)
library(rstan)
library(bayesplot)

source("../pipeline/load_cases_austria.R")

# Settings
prediction_days <- 45
options(mc.cores = 6)

message("Loading movement data")
apple_movement <- tidycovid19::download_apple_mtr_data(type="country_city", cached=TRUE)
apple_movement_austria <- as.data.table(apple_movement)[iso3c=="AUT"][city=="Vienna"]

message("Loading nowcast R estimates")
r_nowcast <- readRDS("../nowcast/austria/W/latest/bigr_estimates.rds")
r_mean <- xts(r_nowcast$mean, r_nowcast$date)

message("Processing mobility data")
driving_raw <- xts(apple_movement_austria$driving/100, ymd(apple_movement_austria$date))
driving_imputed <- na.aggregate(driving_raw, FUN = mean)
driving_raw_l3 <- xts(lag(apple_movement_austria$driving/100, 3), ymd(apple_movement_austria$date))
driving_l3 <- na.aggregate(driving_raw_l3, FUN = mean)
driving_raw_l7 <- xts(lag(apple_movement_austria$driving/100, 7), ymd(apple_movement_austria$date))
driving_l7 <- na.aggregate(driving_raw_l7, FUN = mean)

walking_raw <- xts(apple_movement_austria$walking/100, ymd(apple_movement_austria$date))
walking_imputed <- na.aggregate(walking_raw, FUN = mean)
walking_raw_l3 <- xts(lag(apple_movement_austria$walking/100, 3), ymd(apple_movement_austria$date))
walking_l3 <- na.aggregate(walking_raw_l3, FUN = mean)
walking_raw_l5 <- xts(lag(apple_movement_austria$walking/100, 5), ymd(apple_movement_austria$date))
walking_l5 <- na.aggregate(walking_raw_l5, FUN = mean)
walking_raw_l7 <- xts(lag(apple_movement_austria$walking/100, 7), ymd(apple_movement_austria$date))
walking_l7 <- na.aggregate(walking_raw_l7, FUN = mean)

r0_before <- xts(zooreg(rep(2.3, 30), end = head(time(r_mean) - days(1), 1)))
r_extended <- c(r0_before, r_mean)

glmodel <- stan_model(file = "stan/rt_glm.stan")

merged_data <- merge(
  r_mean, 
  walking_imputed, 
  driving_imputed, 
  all=FALSE)
merged_data_df <- as.data.frame(merged_data)

K <- 14

model_input <- list(
  N0 = match(head(time(merged_data),1), time(walking_imputed)),
  N = length(walking_imputed),
  walking = as.vector(walking_imputed),
  driving = as.vector(driving_imputed),
  
  M = length(merged_data_df$r_mean),
  Rt = merged_data_df$r_mean,
  K = K
)

glm_fit <- sampling(glmodel, 
         data=model_input,
         iter=10000,
         chains=4, 
         control = list(adapt_delta=0.9))

shinystan::launch_shinystan(glm_fit)

extracted_results <- extract(glm_fit,permuted=TRUE)
glm_summary <- summary(glm_fit, probs=c(0.1,0.9))

new_y <- extract(glm_fit, pars="Rt_rep")
pred<-apply(new_y[[1]],2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975)) 
Rt_rep_quantiles <- as.data.frame(t(pred))
Rt_rep_quantiles$date <- time(walking_imputed)[K:length(walking_imputed)]

r_forecast_plot <- ggplot(data=r_nowcast, aes(x=date)) +
  geom_ribbon(aes(ymin=bottom, ymax=top, fill = "Nowcast"), alpha=0.1) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = "Nowcast"), alpha=0.1) + 
  
  geom_ribbon(data=Rt_rep_quantiles, aes(x=ymd(date), ymin=`25%`, ymax=`75%`, fill = "GLM"), alpha=0.1) + 
  geom_ribbon(data=Rt_rep_quantiles, aes(x=ymd(date), ymin=`2.5%`, ymax=`97.5%`, fill = "GLM"), alpha=0.1) + 
  
  geom_hline(yintercept = 1, color="black", size=1, alpha=0.3, linetype="dashed") +
  theme_pubr()

show(r_forecast_plot)
