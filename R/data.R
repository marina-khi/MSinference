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

#' A historical dataset on nominal annual house prices from \cite{Knoll2017}
#' and respective macroeconomic variables from \cite{Jorda2017}.
#' 
#' The dataset consists of nominal house prices and other macroeconomic
#' variables for 8 advanced economies over the time period 1890-2012.
#' 
#' The data are distributed through openICPSR, a public access repository
#' supported by the Inter-university Consortium for
#' Political and Social Research (ICPSR), under an Other License.
#' 
#' The dataset is freely available and was downloaded on 13 January 2022 from
#' https://www.openicpsr.org/openicpsr/project/113055/version/V1/view.
#' 
#' The macroeconomic variables are obtained via
#' the Jordà-Schularick-Taylor Macrohistory Database (\cite{Jorda2017}).
#' The Jordà-Schularick-Taylor Macrohistory Database is freely available and
#' may be used, shared and adapted subject to Creative Commons
#' Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).
#' 
#' The data were downloaded on 13 January 2022 from
#' https://www.macrohistory.net/database/.
#'
#' The entries in the dataset present one of the following measurements:
#' consumer price index (\verb|cpi|), real GDP (\verb|rgdppc|),
#' population size (\verb|pop|), long-term interest rate (\verb|ltrate|),
#' nominal house prices (\verb|hpnom|), real house prices (\verb|hpreal|)
#' obtained through deflation the nominal house prices with the corresponding
#' consumer price index (CPI), inflation (\verb|infl|) which is measured as
#' change in CPI, logarithms of real house prices, population size and real GDP
#' (\verb|log_hp|, \verb|log_gdp| and \verb|log_pop|, respectively) and
#' the changes of the log house prices, log population size, log real GDP,
#' long-term interest rate and inflation (\verb|delta_log_hp|,
#' \verb|delta_log_gdp|, \verb|delta_log_pop|, \verb|delta_ltrate|,
#' and \verb|delta_infl|, respectively). All of the measurements are recorded
#' per country (columns \verb|iso| and \verb|country|) and per year
#' (column \verb|year|).
#' 
#' @format A matrix with 984 rows and 18 columns. Each column corresponds to
#' either the ISO/name of the country and year or one of the measurements
#' described above recorded for that combination of country-year.
#' 
#' @source \url{https://www.openicpsr.org/openicpsr/project/113055/version/V1/view}
#'         \url{https://www.macrohistory.net/database/}
#' @usage data("house_prices")
"house_prices"