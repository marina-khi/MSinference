% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{house_prices}
\alias{house_prices}
<<<<<<< HEAD
\title{House price indices and macro-financial indicators for advanced economies
per country.}
\format{
A matrix with 2214 rows and 18 variables:
\itemize{
  \item iso: ISO code of the country
  \item year: year of the observation
  \item country: name of the country
  \item cpi: consumer price index (CPI) from the Jordà-Schularick-Taylor Macrohistory Database
  \item rgdppc: real GDP from the Jordà-Schularick-Taylor Macrohistory Database
  \item pop: population size from the Jordà-Schularick-Taylor Macrohistory Database
  \item ltrate: long-term interest rate from the Jordà-Schularick-Taylor Macrohistory Database
  \item hpnom: nominal house prices
  \item hpreal: deflated nominal house prices with the CPI
  \item infl: inflation that is measured as change in CPI
  \item log_hp: logarithm of the real house prices
  \item log_gdp: logarithm of the real GDP
  \item log_pop: logarithm of the population size
  \item delta_log_hp: first difference of the logarithm of the real house prices
  \item delta_log_gdp: first difference of the logarithm of the real GDP
  \item delta_log_pop: first difference of the logarithm of the population size
  \item delta_ltrate: first difference of the long-term interest rate
  \item delta_infl: first difference of the inflation
 }
}
\source{
\url{https://www.openicpsr.org/openicpsr/project/113055/version/V1/view},
\url{https://www.macrohistory.net/database/}
=======
\title{A historical dataset on nominal annual house prices from \cite{Knoll2017}
and respective macroeconomic variables from \cite{Jorda2017}.}
\format{
A matrix with 984 rows and 18 columns. Each column corresponds to
either the ISO/name of the country and year or one of the measurements
described above recorded for that combination of country-year.
}
\source{
\url{https://www.openicpsr.org/openicpsr/project/113055/version/V1/view}
        \url{https://www.macrohistory.net/database/}
>>>>>>> fb00f9b83bb9450f09d90344beebc5f014891393
}
\usage{
data("house_prices")
}
\description{
<<<<<<< HEAD
Data that combines residential house price indices  for 14 advanced economies
over the years 1870 to 2012 from Knoll et al. (2017) with several
macro-financial variables from the Jordà-Schularick-Taylor Macrohistory Database.
}
\details{
We .

We also include the logarithms of real house prices (log_hp), the real GDP (log_gdp)
and the population size (log_pop) to measure the growth rates of these variables.

Furthermore, we include the first difference of logarithm of real house
prices (delta_log_hp), the logarithm of the real GDP (delta_log_gdp),
the logarithm of the population size (delta_log_pop), long-term interest
rate (delta_ltrate) and inflation (delta_infl) in order to simplify
further calculations.
=======
The dataset consists of nominal house prices and other macroeconomic
variables for 8 advanced economies over the time period 1890-2012.
}
\details{
The data are distributed through openICPSR, a public access repository
supported by the Inter-university Consortium for
Political and Social Research (ICPSR), under an Other License.

The dataset is freely available and was downloaded on 13 January 2022 from
https://www.openicpsr.org/openicpsr/project/113055/version/V1/view.

The macroeconomic variables are obtained via
the Jordà-Schularick-Taylor Macrohistory Database (\cite{Jorda2017}).
The Jordà-Schularick-Taylor Macrohistory Database is freely available and
may be used, shared and adapted subject to Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).

The data were downloaded on 13 January 2022 from
https://www.macrohistory.net/database/.

The entries in the dataset present one of the following measurements:
consumer price index (\verb|cpi|), real GDP (\verb|rgdppc|),
population size (\verb|pop|), long-term interest rate (\verb|ltrate|),
nominal house prices (\verb|hpnom|), real house prices (\verb|hpreal|)
obtained through deflation the nominal house prices with the corresponding
consumer price index (CPI), inflation (\verb|infl|) which is measured as
change in CPI, logarithms of real house prices, population size and real GDP
(\verb|log_hp|, \verb|log_gdp| and \verb|log_pop|, respectively) and
the changes of the log house prices, log population size, log real GDP,
long-term interest rate and inflation (\verb|delta_log_hp|,
\verb|delta_log_gdp|, \verb|delta_log_pop|, \verb|delta_ltrate|,
and \verb|delta_infl|, respectively). All of the measurements are recorded
per country (columns \verb|iso| and \verb|country|) and per year
(column \verb|year|).
>>>>>>> fb00f9b83bb9450f09d90344beebc5f014891393
}
\keyword{datasets}
