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
#' one coutnry, with the name of the country (denoted by three letter) being
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
