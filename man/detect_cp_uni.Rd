% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{detect_cp_uni}
\alias{detect_cp_uni}
\title{Detect Change Points on an univariate time series.}
\usage{
detect_cp_uni(
  data,
  n_iterations,
  q,
  a = 1,
  b = 1,
  c = 0.1,
  prior_var_phi = 0.1,
  par_theta_c = 1,
  par_theta_d = 1,
  print_progress = TRUE,
  user_seed = 1234L
)
}
\arguments{
\item{data}{vector of observations.}

\item{n_iterations}{number of MCMC iteration.}

\item{q}{probability of performing a split at each iterations.}

\item{a, b, c}{parameters of the Normal-Gamma prior for \eqn{\mu} and \eqn{\lambda}.}

\item{prior_var_phi}{parameters for the correlation coefficient in the likelihood.}

\item{par_theta_c, par_theta_d}{parameters of the shifted Gamma prior for \eqn{\theta}.}

\item{print_progress}{If TRUE (default) print the progress bar.}

\item{user_seed}{seed for random distribution generation.}
}
\value{
Function \code{detect_cp_uni} returns a list containing the following components: \itemize{
\item{\code{$orders}} a matrix where each row corresponds to the output order of the corresponding iteration.
\item{\code{time}} computational time in seconds.
\item{\code{$sigma_MCMC}} traceplot for \eqn{\sigma}.
\item{\code{$sigma_MCMC_01}} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\sigma} was accepted, \eqn{0} otherwise.
\item{\code{$theta_MCMC}} traceplot for \eqn{\theta}.
}
}
\description{
Detect Change Points on an univariate time series.
}
\examples{

data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))

out <- detect_cp_uni(data = data_vec,
                            n_iterations = 2500,
                            q = 0.25)


}
