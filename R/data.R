#' Simulated moult scores
#'
#' A dataset containing simulated moult observations
#'
#' @format A data frame with 500 rows and 3 variables:
#' \describe{
#'   \item{yday}{ordinal_day, such that all active moult observations are continous in the middle of the year}
#'   \item{yday_shifted}{shifted ordinal day, such that the active moult period straddles "new year's eve"}
#'   \item{moult_score}{linear moult progression score [0,1]}
#'   ...
#' }
#'
"sim_data_small"