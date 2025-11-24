#' EU Inflation Dataset (Standardized Matrix Form)
#'
#' This dataset contains standardized monthly inflation rates
#' (annual average rate of change) in the Harmonized Index of Consumer Prices
#' (HICP) for the European Union. Data cover February 1997 to December 2024
#' and are organized across the 12 COICOP expenditure categories.
#'
#' Each row corresponds to one COICOP category, and each column corresponds to
#' a monthly observation.
#'
#' @format A numeric matrix with 12 rows and 355 columns:
#' \describe{
#'   \item{Rows}{COICOP categories}
#'   \item{Columns}{Monthly inflation observations (Feb 1997–Dec 2024)}
#' }
#'
#' @usage data(eu_inflation)
#'
#' @source Eurostat — \emph{HICP (prc_hicp_aind)} available at
#'   the Eurostat Data Explorer.
"eu_inflation"
