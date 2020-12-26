# Needed for Rscript
library(xts)  # to manipulate time series of stock data
library(quantmod)  # to download stock data
library(PerformanceAnalytics)  # to compute performance measures
library(dplyr)
library(CVXR)

# Chunk setting
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  cache = TRUE, 
  fig.align="center",
  fig.pos="t",
  out.width='80%'
)
set.seed(0210)


# options(digits = 3)
# options(dplyr.print_min = 4, dplyr.print_max = 4)


