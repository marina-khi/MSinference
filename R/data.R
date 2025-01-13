#' Number of daily new cases of infections of COVID-19 per country.
#'
#' Data on the geographic distribution of COVID-19 cases worldwide
#' (© ECDC [2005-2019])
#'
#' Each entry in the dataset denotes the number of new cases of infection
#' per day and per country. In order to make the data comparable across
#' countries, we take the day of the 100th confirmed case in each country as
#' the starting date t = 1. This way of “normalizing” the data is
#' common practice (Cohen and Kupferschmidt (2020)).
#' 
#' @format A matrix with 99 rows and 41 columns. Each column corresponds to
#' one country, with the name of the country (denoted by three letter) being
#' the name of the column.
#' @source \url{https://www.ecdc.europa.eu/en}
#' @usage data("covid")
"covid"

#' House price indices and macro-financial indicators for advanced economies
#' per country.
#'
#' Data that combines residential house price indices  for 14 advanced economies
#' over the years 1870 to 2012 from Knoll et al. (2017) with several
#' macro-financial variables from the Jordà-Schularick-Taylor Macrohistory Database.
#'
#' We .
#' 
#' We also include the logarithms of real house prices (log_hp), the real GDP (log_gdp)
#' and the population size (log_pop) to measure the growth rates of these variables.
#' 
#' Furthermore, we include the first difference of logarithm of real house
#' prices (delta_log_hp), the logarithm of the real GDP (delta_log_gdp),
#' the logarithm of the population size (delta_log_pop), long-term interest
#' rate (delta_ltrate) and inflation (delta_infl) in order to simplify
#' further calculations.
#' 
#' @format A matrix with 2214 rows and 18 variables:
#' \itemize{
#'   \item iso: ISO code of the country
#'   \item year: year of the observation
#'   \item country: name of the country
#'   \item cpi: consumer price index (CPI) from the Jordà-Schularick-Taylor Macrohistory Database
#'   \item rgdppc: real GDP from the Jordà-Schularick-Taylor Macrohistory Database
#'   \item pop: population size from the Jordà-Schularick-Taylor Macrohistory Database
#'   \item ltrate: long-term interest rate from the Jordà-Schularick-Taylor Macrohistory Database
#'   \item hpnom: nominal house prices
#'   \item hpreal: deflated nominal house prices with the CPI
#'   \item infl: inflation that is measured as change in CPI
#'   \item log_hp: logarithm of the real house prices
#'   \item log_gdp: logarithm of the real GDP
#'   \item log_pop: logarithm of the population size
#'   \item delta_log_hp: first difference of the logarithm of the real house prices
#'   \item delta_log_gdp: first difference of the logarithm of the real GDP
#'   \item delta_log_pop: first difference of the logarithm of the population size
#'   \item delta_ltrate: first difference of the long-term interest rate
#'   \item delta_infl: first difference of the inflation
#'  }
#' @source \url{https://www.openicpsr.org/openicpsr/project/113055/version/V1/view},
#' \url{https://www.macrohistory.net/database/}
#' @usage data("house_prices")
"house_prices"

#' Hadley Centre Central England Temperature (HadCET) dataset,
#' Monthly Mean Central England Temperature (Degrees C)
#'
#' The CET dataset is the longest instrumental record of temperature
#' in the world. It contains the mean monthly surface air temperatures
#' (in degrees Celsius) from the year 1659 to the present. These monthly
#' temperatures are representative of a roughly triangular area of
#' the United Kingdom enclosed by Lancashire, London and Bristol.
#' Manley (1953, 1974) compiled most of the monthly series,
#' covering 1659 to 1973.  These data were updated to 1991 by
#' Parker et al (1992). It is now kept up to date by
#' the Climate Data Monitoring section of the Hadley Centre, Met Office.
#' 
#' Since 1974 the data have been adjusted to allow for urban warming:
#' currently a correction of -0.2 C is applied to mean temperatures.
#' CET datasets are freely available for use under Open Government License.
#' 
#' @format A numeric vector of length 359.
#' @source \url{https://www.metoffice.gov.uk/hadobs/hadcet/}
#' @usage data("temperature")
"temperature"
