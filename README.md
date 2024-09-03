# MSinference: An R Package for the Multiscale Inference for Nonparametric Time Trend(s)

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/MSinference)](https://CRAN.R-project.org/package=MSinference)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
[![](https://cranlogs.r-pkg.org/badges/grand-total/MSinference)](https://CRAN.R-project.org/package=MSinference)
<!-- badges: end -->

Author: Marina Khismatullina

## Description
The package performs a multiscale analysis of a nonparametric regression or nonparametric regressions with time series errors. In case of one regression, with the help of this package it is possible to detect the regions where the trend function is increasing or decreasing. In case of multiple regressions, the test identifies regions where the trend functions are different from each other. See Khismatullina and Vogt (2020) <doi:10.1111/rssb.12347>, Khismatullina and Vogt (2022) <doi:10.48550/arXiv.2209.10841> and Khismatullina and Vogt (2023) <doi:10.1016/j.jeconom.2021.04.010> for more details on theory and applications.

## Dependencies

To optimize computational performance, most methods implemented in the package are written in C++. However, for user convenience, the interface of the package is entirely implemented in R, with minimal dependencies, including:

- `Rcpp`: facilitates seamless integration between R and C++.
- `foreach`, `parallel`, `doParallel`: are necessary for fast parallel computation of the critical values.

## Issues, bug reports, contributions, further help

You can raise issues, report bugs, seek for further help, or submit your contribution to the R package `MSinference` at the
github repository [marina-khi/MSinference](https://github.com/marina-khi/MSinference). 

For bug reports, you are kindly asked to make a small and self-contained program which exposes the bug.