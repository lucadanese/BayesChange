#' Multivariate Synthetic Epidemiological Time-Series Data
#'
#' A multivariate synthetic infection-count dataset generated from three
#' independent stochastic epidemic processes simulated using the
#' Doob–Gillespie algorithm. Each time series has a different
#' transmission rate, resulting in distinct change point
#' structures.
#'
#' The simulation follows the stochastic framework described in:
#' Anderson, D. F. and Kurtz, T. G. (2015).
#' *Stochastic Analysis of Biochemical Systems*.
#' Springer International Publishing.
#'
#' @format A \eqn{3 \times 200} numeric matrix.
#' Each row represents one synthetic epidemic time series and each column
#' corresponds to a discrete time point.
#'
#' @details
#' All three epidemic processes use:
#' \itemize{
#'   \item \code{S0 = 100000}, \code{I0 = 20}
#'   \item \code{max_time = 200}
#'   \item \code{xi_0 = 1/8}
#' }
#'
#' The three \code{beta} vectors differ in their change-point locations and
#' values:
#' \itemize{
#'   \item Series 1: \code{0.211 → 0.55} at time 120
#'   \item Series 2: \code{0.215 → 0.52} at time 120
#'   \item Series 3: \code{0.193 → 0.53} at time 30
#' }
#'
#' @usage data(epi_synthetic_multi)
#'
"epi_synthetic_multi"
