# External Imports
require(data.table, quietly = TRUE) 
require(future, quietly = TRUE)
require(lubridate)

# Internal imports
source("../pipeline/load_cases_austria.R")

## Require for data and nowcasting
# require(EpiNow, quietly = TRUE)
# require(NCoVUtils, quietly = TRUE)

## Required for forecasting
# require(future.apply, quietly = TRUE)
# require(fable, quietly = TRUE)
# require(fabletools, quietly = TRUE)
# require(feasts, quietly = TRUE)
# require(urca, quietly = TRUE)
NCoVUtils::reset_cache()

# Load delays
delay_defs <- readRDS("../delays/data/delay_defs.rds")
incubation_defs <- readRDS("../delays/data/incubation_defs.rds")

# Load cases
cases <- load_cases_austria()

# Load cases for Austria ------------------------------------------
message("Loading regional case data for Austria")


# Run pipeline ----------------------------------------------------
future::plan("multiprocess", workers = future::availableCores())

message("Running regional rt pipeline")
EpiNow::regional_rt_pipeline(
  cases = cases,
  delay_defs = delay_defs,
  incubation_defs = incubation_defs,
  target_folder = "austria",
  case_limit = 60,
  horizon = 14,
  nowcast_lag = 8,
  rt_samples = 15,
  approx_delay = TRUE,
  report_forecast = TRUE, 
  verbose = TRUE,
  forecast_model = function(...) {
    EpiSoon::fable_model(model = fabletools::combination_model(fable::RW(y ~ drift()), fable::ETS(y), 
                                                               fable::NAIVE(y),
                                                               cmbn_args = list(weights = "inv_var")), ...)
  }
)
