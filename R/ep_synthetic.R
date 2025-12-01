#' Synthetic Epidemiological Time-Series Data
#'
#' A toy dataset generated from a stochastic SIR-type epidemic model using
#' the Doob–Gillespie algorithm. The transmission rate \code{beta} changes
#' once over time, resultiing in a single change-point in the infection
#' counts.
#'
#' The simulation follows the stochastic simulation framework of:
#' Anderson, D. F. and Kurtz, T. G. (2015). *Stochastic Analysis of
#' Biochemical Systems*. Springer International Publishing.
#'
#' @format A 200 × 1 numeric matrix containing the daily number of infection events.
#'
#' @details
#' The simulation uses:
#' \itemize{
#'   \item \code{S0 = 10000}, \code{I0 = 50}
#'   \item \code{max_time = 200}
#'   \item A piecewise-constant transmission rate vector with a change at time 130
#'   \item Infection event times aggregated using \code{floor()}
#' }
#'
#' @usage data(epi_synthetic)
#'
"epi_synthetic"
