#' Daily Stock Price Dataset (Multivariate Open/Close Series)
#'
#' This dataset contains detrended and standardized daily opening and closing
#' stock prices for the 50 largest companies in the S&P 500 index. Data cover
#' the period from January 1, 2020 to January 1, 2022, with 505 daily
#' observations per company.
#'
#' @format A numeric array of dimension \code{2 × 505 × 50}:
#' \describe{
#'   \item{Dimension 1}{Price type: \code{"Open"}, \code{"Close"}}
#'   \item{Dimension 2}{Daily observations (505)}
#'   \item{Dimension 3}{Companies (S&P 500 top 50)}
#' }
#'
#' @usage data(stock_multi)
#'
#' @source Yahoo Finance (https://finance.yahoo.com/)
"stock_multi"
