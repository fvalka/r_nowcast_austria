# Packages -----------------------------------------------------------------
require(data.table, quietly = TRUE) 
require(future, quietly = TRUE)
require(lubridate)
## Require for data and nowcasting
# require(EpiNow, quietly = TRUE)
# require(NCoVUtils, quietly = TRUE)

## Required for forecasting
# require(future.apply, quietly = TRUE)
# require(fable, quietly = TRUE)
# require(fabletools, quietly = TRUE)
# require(feasts, quietly = TRUE)
# require(urca, quietly = TRUE)

# Get cases ---------------------------------------------------------------

NCoVUtils::reset_cache()

# Get linelist ------------------------------------------------------------

message("Loading line-list")
linelist <- 
  data.table::fread("https://raw.githubusercontent.com/epiforecasts/NCoVUtils/master/data-raw/linelist.csv")

#linelist <- get_linelist()


delays <- linelist[!is.na(date_onset_symptoms)][, 
                                                .(report_delay = as.numeric(lubridate::dmy(date_confirmation) - 
                                                                              as.Date(lubridate::dmy(date_onset_symptoms))))]

delays <- delays$report_delay

# Set up cores -----------------------------------------------------
if (!interactive()){
  options(future.fork.enable = TRUE)
}


future::plan("multiprocess", workers = round(future::availableCores() / 3))

# Fit the reporting delay -------------------------------------------------

message("Calculating delay distribution")
delay_defs <- EpiNow::get_dist_def(delays,
                                   bootstraps = 100, 
                                   samples = 1000)

saveRDS(delay_defs, "data/delay_defs.rds")


## Mean delay
mean_delay <- exp(EpiNow::covid_incubation_period[1, ]$mean)
message(sprintf("Mean delay: %f days",mean_delay))