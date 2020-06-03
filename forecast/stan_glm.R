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
plan(multiprocess)

source("../pipeline/load_cases_austria.R")

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

message("Run STAN GLM fit")
fit <- stan_glm(r_mean ~ 
                  walking_imputed + walking_l5 + walking_l7 +
                  #driving_imputed + driving_l3 + 
                  driving_l7,
                data=as.data.frame(merge(
                  r_mean, 
                  walking_imputed, walking_l3, walking_l5, walking_l7, 
                  driving_imputed, driving_l3, driving_l7,
                  all=FALSE)))

plot(fit, prob = 0.5)
prior_summary(fit)

message("Performing posterior prediction on mobility data")
pprediction <- posterior_predict(fit, 
                                 as.data.frame(merge(
                                   walking_imputed, walking_l3, walking_l5, walking_l7, 
                                   driving_imputed, driving_l3, driving_l7,
                                   all=FALSE)))

posterior_quantiles <- gather(as.data.frame(pprediction), "date", "Rt") %>%
  dplyr::group_by(date) %>% summarise(
    cri95_lower=quantile(Rt, 0.025),
    cri95_upper=quantile(Rt, 0.975),
    cri50_lower=quantile(Rt, 0.25),
    cri50_upper=quantile(Rt, 0.75),
    Rt_mean = mean(Rt))

r_forecast_plot <- ggplot(data=r_nowcast, aes(x=date)) +
  geom_ribbon(aes(ymin=bottom, ymax=top, fill = "Nowcast"), alpha=0.1) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = "Nowcast"), alpha=0.1) + 
  
  geom_ribbon(data=posterior_quantiles, aes(x=ymd(date), ymin=cri50_lower, ymax=cri50_upper, fill = "GLM"), alpha=0.1) + 
  geom_ribbon(data=posterior_quantiles, aes(x=ymd(date), ymin=cri95_lower, ymax=cri95_upper, fill = "GLM"), alpha=0.1) + 
  
  geom_hline(yintercept = 1, color="black", size=1, alpha=0.3, linetype="dashed") +
  theme_pubr()

ggsave("figures/latest/r-forecast.pdf", plot=r_forecast_plot, width=16, height=9)


###############################
# GLM Case Forecast Simulation
###############################

message("Performing case prediction based upon R prediction")
simSize <- 1000

message("Loading and processing case nowcast")
nowcast_data <- readRDS("../nowcast/austria/W/latest/nowcast.rds")

nowcast_mean <- nowcast_data %>% group_by(date) %>% summarise(
  cri95_lower=quantile(cases, 0.025),
  cri95_upper=quantile(cases, 0.975),
  cri50_lower=quantile(cases, 0.25),
  cri50_upper=quantile(cases, 0.75),
  cases = mean(cases))

message("Extend R prediction")
dates_ppred <- ymd(names(pprediction[1,]))
ppred_extended <- future_map(1:dim(pprediction)[1], ~ c(zoo(pprediction[.x,], dates_ppred), 
                        zooreg(rep(mean(tail(pprediction[.x,],7)), 14), start = tail(dates_ppred,1) + days(1))))


message("Predicting cases using monte carlo simulation based on nowcast, R prediction and generation intervals")
tic()
cases_pred_r_nowcast_raw <- future_map_dfr(1:simSize, ~ predict_cases(cases = nowcast_data[sample==sample(1:1000,1)],
                                                                    rts = data.frame("date"=dates_ppred, "rt"=pprediction[sample(1:dim(pprediction)[1],1),]),
                                                                    forecast_date = as.Date("2020-03-1"),
                                                                    serial_interval = EpiNow::covid_generation_times[,sample(1:dim(EpiNow::covid_generation_times)[2],1)]),
                                           .progress = TRUE)
toc()

cases_pred_r_nowcast_raw$sample <- ceiling(index(cases_pred_r_nowcast_raw)/length(unique(cases_pred_r_nowcast_raw$date)))

cases_ppred <- cases_pred_r_nowcast_raw %>%
  dplyr::group_by(date) %>% summarise(
    cri95_lower=quantile(cases, 0.025),
    cri95_upper=quantile(cases, 0.975),
    cri50_lower=quantile(cases, 0.25),
    cri50_upper=quantile(cases, 0.75),
    cases = mean(cases))

cases_austria <- load_cases_austria()
local_cases_vienna <- cases_austria[region=="W"] %>% group_by(date) %>% summarise(cases=sum(confirm))

cases_forecast_plot <- ggplot(data=cases_ppred, aes(x=date)) +
  geom_ribbon(data=nowcast_mean, aes(x=date, ymin=cri95_lower, ymax=cri95_upper, fill="nowcast"), alpha=0.2) +
  geom_ribbon(data=nowcast_mean, aes(x=date, ymin=cri50_lower, ymax=cri50_upper, fill="nowcast"), alpha=0.2) +
  
  geom_ribbon(data=cases_ppred, aes(x=date, ymin=cri95_lower, ymax=cri95_upper, fill="cases_ppred"), alpha=0.2) +
  geom_ribbon(data=cases_ppred, aes(x=date, ymin=cri50_lower, ymax=cri50_upper, fill="cases_ppred"), alpha=0.2) +
  
  geom_point(data=local_cases_vienna, aes(x=date, y=cases)) +
  
  ggtitle("Case Prediction") +
  #coord_cartesian(ylim = c(0, 200)) +
  theme_pubr()


ggsave("figures/latest/cases-forecast.pdf", plot=cases_forecast_plot, width=16, height=9)


message("Simulating incubation time and reporting delay on predicted cases")
# Load delays
delay_defs <- readRDS("../delays/data/delay_defs.rds")
incubation_defs <- readRDS("../delays/data/incubation_defs.rds")

delay_def <- EpiNow::lognorm_dist_def(mean = 5, mean_sd = 1,
                                      sd = 3, sd_sd = 1, max_value = 30,
                                      samples = 1, to_log = TRUE)

cases_ppred_rounded <- data.frame("date"=cases_ppred$date, "cases"=round(cases_ppred$cases))

# Simulate Delay
tic()
simulated_reported_cases <- future_map_dfr(1:10000, ~ EpiNow::adjust_infection_to_report(
  cases_pred_r_nowcast_raw[sample==sample(1:simSize,1)],
  delay_def = delay_defs[sample(1:dim(delay_defs)[1],1)],
  incubation_def = incubation_defs[sample(1:dim(incubation_defs)[1],1)], 
  reporting_effect = rep(1, 7),
  type = "sample", return_onset = FALSE),
  .progress = TRUE)
toc()

simulated_reported_cases_stats <- simulated_reported_cases %>%
  dplyr::group_by(date) %>% summarise(
    cri95_lower=quantile(cases, 0.025),
    cri95_upper=quantile(cases, 0.975),
    cri50_lower=quantile(cases, 0.25),
    cri50_upper=quantile(cases, 0.75),
    cases = mean(cases))

cases_delayed_plot <- ggplot(data=cases_ppred, aes(x=date)) +
  geom_ribbon(data=nowcast_mean, aes(x=date, ymin=cri95_lower, ymax=cri95_upper, fill="nowcast"), alpha=0.2) +
  geom_ribbon(data=nowcast_mean, aes(x=date, ymin=cri50_lower, ymax=cri50_upper, fill="nowcast"), alpha=0.2) +
  
  geom_ribbon(data=simulated_reported_cases_stats, aes(x=date, ymin=cri95_lower, ymax=cri95_upper, fill="cases_ppred"), alpha=0.2) +
  geom_ribbon(data=simulated_reported_cases_stats, aes(x=date, ymin=cri50_lower, ymax=cri50_upper, fill="cases_ppred"), alpha=0.2) +
  
  geom_point(data=local_cases_vienna, aes(x=date, y=cases)) +
  
  ggtitle("Case Prediction") +
  theme_pubr()


ggsave("figures/latest/cases-forecast-delayed.pdf", plot=cases_delayed_plot, width=16, height=9)

