library(tidyverse)
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
plan(multisession)

source("../pipeline/load_cases_austria.R")

message("Simulating incubation time and reporting delay on predicted cases")
# Load delays
delay_defs <- readRDS("../delays/data/delay_defs.rds")
incubation_defs <- readRDS("../delays/data/incubation_defs.rds")

# Load simulated cases 
cases_pred_r_nowcast_raw <- readRDS("output/raw-cases-prediction.rds")
cases_pred_r_nowcast_raw$sample <- ceiling(index(cases_pred_r_nowcast_raw)/length(unique(cases_pred_r_nowcast_raw$date)))

# Load actual cases
cases_austria <- load_cases_austria("../nowcast/data/cases.csv")
local_cases_vienna <- cases_austria[region=="W"] %>% group_by(date) %>% summarise(cases=sum(confirm))

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
saveRDS(simulated_reported_cases, "output/simulated_reported_cases.rds")

simulated_reported_cases_stats <- simulated_reported_cases %>%
  dplyr::group_by(date) %>% summarise(
    cri95_lower=quantile(cases, 0.025),
    cri95_upper=quantile(cases, 0.975),
    cri50_lower=quantile(cases, 0.25),
    cri50_upper=quantile(cases, 0.75),
    cases = mean(cases))

colors <- c("Prediction"="#70c1b3", "Actual Cases"="#247ba0")

cases_delayed_plot <- ggplot(data=simulated_reported_cases_stats, aes(x=date)) +
  
  geom_ribbon(aes(x=date, ymin=cri95_lower, ymax=cri95_upper, fill="Prediction"), alpha=0.2) +
  geom_ribbon(aes(x=date, ymin=cri50_lower, ymax=cri50_upper, fill="Prediction"), alpha=0.2) +
  
  geom_point(data=local_cases_vienna, aes(x=date, y=cases, color="Actual Cases")) +
  #geom_ma(data=local_cases_vienna, aes(x=date - days(3), y=cases, color="Actual Cases"), ma_fun = SMA, n=7) + 
  
  ggtitle("Case Prediction") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(color="", fill="") + 
  theme_pubr()


ggsave("figures/latest/cases-forecast-delayed.pdf", plot=cases_delayed_plot, width=16, height=9)
