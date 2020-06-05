require(lubridate, quietly = TRUE)
require(data.table, quietly = TRUE)

load_cases_austria <- function(file = "https://raw.githubusercontent.com/osaukh/dashcoch-AT/master/data_AT/covid19_cases_austria.csv") {
  d_Austria_per_county <- read.csv(file, 
                                   fileEncoding="UTF-8-BOM")
  
  dates <- ymd(d_Austria_per_county$Date[-1])
  cases <- data.table(date=dates)
  
  states <- c("W", "NÖ", "OÖ", "T", "K", "ST", "S", "V", "B", "AT")
  #states <- c("W")
  
  for (state in states) {
    # CSH/Github pro Bundesland
    cases_state <- diff(d_Austria_per_county[,state])
    
    # Some states have negative cases, these will be applied to the previous day and set to 0
    repeat{
      found_negative_cases <- FALSE
      for (i in 2:length(cases_state)) {
        if(cases_state[i] < 0){
          print(sprintf("Warning negative cases found on day %s correcting by applying those corrections to the previous day", dates[i]))
          cases_state[i-1] <- cases_state[i-1] + cases_state[i]
          cases_state[i] <- 0
          found_negative_cases <- TRUE
        }
      }
      if (!found_negative_cases) { break }
    }
    
    cases[, (state) := cases_state]
  }
  
  cases <- data.table::melt(as.data.table(cases), id.vars=c("date"), variable.name="region", value.name="cases")
  cases$date <- ymd(cases$date)
  
  cases <- data.table::setDT(cases)[!is.na(region)][, 
                                                    `:=`(local = cases, imported = 0)][, cases := NULL]
  
  cases <- data.table::melt(cases, measure.vars = c("local", "imported"),
                            variable.name = "import_status",
                            value.name = "confirm")
  
  return(cases)
}