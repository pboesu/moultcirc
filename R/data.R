#' Simulated moult scores
#'
#' A dataset containing simulated moult observations
#'
#' @format A data frame with 500 rows and 5 variables:
#' \describe{
#'   \item{yday}{ordinal_day, such that all active moult observations are continous in the middle of the year}
#'   \item{yday_shifted}{shifted ordinal day, such that the active moult period straddles "new year's eve"}
#'   \item{moult_score}{linear moult progression score [0,1]}
#'   \item{moult_cat}{moult category, Factor with 3 levels: "O","M","N"}
#'   \item{moult_numcat}{moult category, numeric: 1,2,3}
#'   ...
#' }
#'
"sim_data_small"
#' Simulated moult scores
#'
#' A dataset containing simulated moult observations
#'
#' @format A data frame with 5000 rows and 5 variables:
#' \describe{
#'   \item{yday}{ordinal_day, such that all active moult observations are continous in the middle of the year}
#'   \item{yday_shifted}{shifted ordinal day, such that the active moult period straddles "new year's eve"}
#'   \item{moult_score}{linear moult progression score [0,1]}
#'   \item{moult_cat}{moult category, Factor with 3 levels: "O","M","N"}
#'   \item{moult_numcat}{moult category, numeric: 1,2,3}
#'   ...
#' }
#'
"sim_data"
