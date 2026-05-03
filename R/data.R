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

#' Major League Baseball No-Hitters
#'
#' A dataset containing a list of Major League Baseball no-hitter games,
#' compiled from Retrosheet, Inc. (\url{http://www.retrosheet.org}).
#'
#' @format A data frame with \code{n} observations and 6 variables:
#' \describe{
#'   \item{Date}{Character. The date of the game (MM/DD/YYYY format).}
#'   \item{Pitcher}{Character. Name of the pitcher who threw the no-hitter.}
#'   \item{Team}{Character. Team for which the pitcher played.}
#'   \item{Lg}{Character. League of the team (e.g., NA, NL, AL). May contain missing values.}
#'   \item{Plate.Umpire}{Character. Name of the home plate umpire.}
#'   \item{Notes}{Character. Additional notes (e.g., "Perfect Gm"). May be empty.}
#' }
#'
#' @details
#' The \code{no_hitters} dataset provides historical records of no-hitter
#' games in Major League Baseball. A no-hitter is a game in which a pitcher
#' (or multiple pitchers) allows no hits over the course of a complete game.
#'
#' The data are sourced from Retrosheet and include early baseball history,
#' with records dating back to the 19th century.
#'
#' @source
#' Retrosheet, Inc. \url{http://www.retrosheet.org}
#'
#' @examples
#' data(no_hitters)
#' head(no_hitters)
#'
"no_hitters"
