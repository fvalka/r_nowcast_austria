# Nowcast of the Time-Varying Reproduction Number for Austria

This nowcast applies the [EpiNow](https://github.com/epiforecasts/EpiNow) package 
to regional COVID-19 case data for Austria. 

***Warning: Work in progress!***

## Data Source
### Austrian Regional Case Data
Austrian regional data is taken from the 
Complexity Science Hubs
[dashboard](https://github.com/osaukh/dashcoch-AT/). They provide a time series of 
daily regional updates from the 
[Bundesministerium für Soziales, Gesundheit, Pflege und Konsumentenschutz](https://www.sozialministerium.at/Informationen-zum-Coronavirus/Neuartiges-Coronavirus-(2019-nCov).html).

### Time Delays
Incubation times are taking from the 
[EpiNow](https://github.com/epiforecasts/EpiNow) package and line-listings cachend in 
[NCovUtils](https://github.com/epiforecasts/NCoVUtils) are currently being used. 
