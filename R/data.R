#' Weekly Beverage Sales in the U.S.
#'
#'
#' @format A data frame with 391 rows and 19 variables:
#' \describe{
#'   \item{STORE}{Store identifier.}
#'   \item{UPC}{Product UPC code.}
#'   \item{date}{Observation date.}
#'   \item{MOVE}{Observed sales count.}
#'   \item{Buy}{Purchase indicator.}
#' }
#'
#' @details
#' Weekly count series of product sales at Dominick's Finer Foods, a now defunct
#' U.S. grocery chain that operated in Chicago, IL and adjacent
#' areas from 1918 to 2013 (universal product code 4640055081, store 81).
#'
#' @source https://www.chicagobooth.edu/research/kilts/research-data/dominicks.
#'
#' @examples
#' data(drinksales)
#' str(drinksales)
#' drinksales$date <- as.Date(drinksales$date, format =  "%d%b%Y")
#' plot(drinksales$date, drinksales$MOVE, type = "l")
"drinksales"
