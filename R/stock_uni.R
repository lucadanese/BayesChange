#' Daily Stock Price Dataset (Univariate Mean Series)
#'
#' This dataset contains detrended and standardized daily mean stock prices for
#' the 50 largest companies (by market capitalization) in the S&P 500 index.
#' Data cover the period from January 1, 2020 to January 1, 2022.
#'
#' For each company, the mean price is computed as the average between the
#' daily high and low prices. Each resulting time series contains 505 daily
#' observations.
#'
#' @format A numeric matrix with 50 rows and 505 columns:
#' \describe{
#'   \item{Rows}{Companies (S&P 500 top 50)}
#'   \item{Columns}{Daily mean prices (505 observations)}
#' }
#'
#' @usage data(stock_uni)
#'
#' @source Yahoo Finance (https://finance.yahoo.com/)
"stock_uni"
