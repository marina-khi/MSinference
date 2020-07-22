## code to prepare `covid` dataset goes here

library(tidyr)
covid_tmp         <- read.csv("/Users/missius/Desktop/Work/multiscale_inference/Epidemic/Code/data/covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "")
covid_tmp         <- covid_tmp[!is.na(covid_tmp$countryterritoryCode), ]
covid_tmp$dateRep <- as.Date(covid_tmp$dateRep, format = "%d/%m/%Y")
covid_tmp         <- complete(covid_tmp, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'),
                              countryterritoryCode, fill = list(cases = 0, deaths = 0))

covid_tmp$cumcases        <- 0
covid_tmp$cumdeaths       <- 0

covid_list <- list()
for (country in unique(covid_tmp$countryterritoryCode)){
  covid_tmp[covid_tmp$countryterritoryCode == country, "cumcases"] <- cumsum(covid_tmp[covid_tmp$countryterritoryCode == country, "cases"])
  covid_tmp[covid_tmp$countryterritoryCode == country, "cumdeaths"] <- cumsum(covid_tmp[covid_tmp$countryterritoryCode == country, "deaths"])
  tmp <- max(covid_tmp[covid_tmp$countryterritoryCode == country, "cumdeaths"])
  if (tmp >= 1000){
    #We restrict our attention only to the contries with more than 1000 deaths and only starting from 100th case
    covid_list[[country]] <- covid_tmp[(covid_tmp$countryterritoryCode == country & covid_tmp$cumcases >= 100),
                                       c("dateRep", "cases", "deaths", "cumcases", "cumdeaths")]
  }
}

max_dates <- max(sapply(covid_list[names(covid_list) != "CHN"], NROW))

countries     <- names(covid_list)
dates         <- unique(covid_tmp$dateRep)
countries_num <- length(covid_list)

covid           <- matrix(NA, ncol = countries_num, nrow = max_dates)
colnames(covid) <- countries

i = 1
for (country in countries) {
  covid[, i] <- covid_list[[country]]$cases[1:max_dates]
  i = i + 1
}

usethis::use_data(covid, overwrite = TRUE)

## code to prepare `temperature` dataset goes here
temperature <- read.table("/Users/missius/Desktop/Work/multiscale_inference/One_trend/Code/data/cetml1659on.dat", header = TRUE, skip = 6)
temperature <- temperature[temperature$YEAR > -99, 'YEAR']
usethis::use_data(temperature, overwrite = TRUE)
