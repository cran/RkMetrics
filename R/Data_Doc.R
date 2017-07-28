#' Switzerland Mortality Data
#'
#' Exposed to Risk and number of deaths data.
#'
#' Mortality data for both Males and Females in Switzerland, from 1981
#' to 2014.
#'
#' These data are freely available at the Human Mortality Database
#'
#' @docType data
#'
#'
#' @keywords datasets
#'
#' @format A data frame with 6 columns corresponding to:
#' \describe{
#'   \item{Year}{Corresponding year of , in US dollars}
#'   \item{Age}{Age of the individual}
#'   \item{E.Male}{Male Exposed-to-Risk Population}
#'   \item{E.Female}{Female Exposed-to-Risk Population}
#'   \item{D.Male}{Number of male death counts, for the given year and age}
#'   \item{D.Female}{Number of female death counts, for the given year and age}
#' }
#'
#' @source \url{http://www.mortality.org/cgi-bin/hmd/country.php?cntr=CHE&level=1}
#'
#' @references Glei, D. and Andreeva, M. (2016). About mortality data for switzerland.
#' 
"Mortality"