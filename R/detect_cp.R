#' @name detect_cp
#' @export detect_cp
#'
#' @title Detect change points on time series.
#' @description The \code{detect_cp} function detect change points on univariate and multivariate time series.
#'
#'
#' @param data if kernel = "ts" a vector or a matrix. If kernel = "epi" a matrix.
#'
#' @param n_iterations number of MCMC iterations.
#' @param n_burnin number of iterations that must be excluded when computing the posterior estimate.
#' @param q probability of performing a split at each iteration.
#' @param kernel can be "ts" if data are time series or "epi" if data are epidemic diffusions.
#' @param print_progress If TRUE (default) print the progress bar.
#' @param user_seed seed for random number generation (default 1234). The output is fully determined by \code{user_seed}: two calls with the same data, arguments and \code{user_seed} return identical results, regardless of any previous call to \code{set.seed()}. The state of R's random number generator is restored on exit. If \code{NULL}, the current state of R's random number generator is used.
#'
#'
#' @param params the model parameters, created with [params_uni()] (univariate
#'   time series), [params_multi()] (multivariate time series) or [params_epi()]
#'   (epidemic diffusions). Parameters not given take their default values, see
#'   [params_uni()]. A plain list with the same names is also accepted.
#'
#' @return A \code{DetectCpObj} object containing:
#'
#' \itemize{
#'   \item \code{$data} Vector or matrix with the input data.
#'   \item \code{$n_iterations} Total number of MCMC iterations.
#'   \item \code{$n_burnin} Number of burn-in iterations removed from posterior summaries.
#'   \item \code{$orders} Matrix of change–point allocations. Each row corresponds to an iteration; columns correspond to time indices.
#'   \item \code{$time} Computational time (in seconds).
#'
#'   \item \code{$entropy_MCMC} A \code{coda::mcmc} object containing the MCMC samples of the entropy measure.
#'   \item \code{$lkl_MCMC} A \code{coda::mcmc} object containing the MCMC samples of the log-likelihood.
#'   \item \code{$phi_MCMC} A \code{coda::mcmc} object storing MCMC draws of \eqn{\gamma}.
#'   \item \code{$sigma_MCMC} A \code{coda::mcmc} object storing MCMC draws of \eqn{\sigma}.
#'   \item \code{$delta_MCMC} A \code{coda::mcmc} object storing MCMC draws of \eqn{\delta}.
#'   \item \code{$I0_MCMC} A \code{coda::mcmc} object storing MCMC draws of \eqn{I_0}.
#'
#'   \item \code{$kernel_ts} Logical; TRUE if the data are time series.
#'   \item \code{$kernel_epi} Logical; TRUE if the data follow an epidemic diffusion process.
#'   \item \code{$univariate_ts} Logical; TRUE if the time series is univariate, FALSE for multivariate.
#' }

#'
#' @examples
#'
#' ## Univariate time series
#'
#' data("eu_inflation")
#'
#' out <- detect_cp(data = eu_inflation[1,], n_iterations = 2000,
#'                  n_burnin = 250, q = 0.5,
#'                  params = params_uni(c = 1, prior_var_phi = 0.1),
#'                  kernel = "ts")
#'
#' print(out)
#'
#' ## Multivariate time series
#'
#' data("eu_inflation")
#'
#' out <- detect_cp(data = eu_inflation[1:3,], n_iterations = 2000,
#'                  n_burnin = 250, q = 0.5,
#'                  params = params_multi(k_0 = 1, nu_0 = 10), kernel = "ts")
#'
#' print(out)
#'
#' \donttest{
#' ## Epidemic diffusions
#'
#' data("epi_synthetic")
#'
#' out <- detect_cp(data = epi_synthetic, n_iterations = 2000, n_burnin = 250,
#'                  q = 0.25, params = params_epi(M = 100, I0_var = 0.1),
#'                  kernel = "epi")
#'
#' print(out)
#' }
#'
#' @references
#'
#' Martínez, A. F., & Mena, R. H. (2014). On a Nonparametric Change Point Detection Model in Markovian Regimes. \emph{Bayesian Analysis}, 9(4), 823–858. \doi{10.1214/14-BA878}
#'
#' Corradin, R., Danese, L., & Ongaro, A. (2022). Bayesian nonparametric change point detection for multivariate time series with missing observations. \emph{International Journal of Approximate Reasoning}, 143, 26--43. \doi{10.1016/j.ijar.2021.12.019}
#'
detect_cp <- function(data,
                      n_iterations,
                      n_burnin = 0,
                      q = 0.5,
                      params = list(),
                      kernel,
                      print_progress = TRUE,
                      user_seed = 1234){

  if((!is.vector(data)) && (!is.matrix(data))) stop("data must be a vector or a matrix")
  if(is.null(n_iterations)) stop("missing number of iterations")
  if(n_iterations < 1) stop("number of iterations must be positive and > 1")
  if(n_burnin < 0) stop("number of burn in iterations must be positive and > 1")
  if((!is.null(q)) && ((q >= 1) | (q <= 0))) stop("q must be in (0,1)")
  if((is.null(kernel) | !(kernel %in% c("ts", "epi")))) stop("kernel must be 'ts' or 'epi'")
  if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
  if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

  user_seed <- local_user_seed(user_seed)

  if(kernel == "ts" && is.vector(data)){

    rp <- resolve_params(params, "uni")
    report_sampler_settings(rp$sampler, "detect_cp")
    pr <- rp$values

    out <- detect_cp_uni(data = data,
                         n_iterations = n_iterations,
                         q = q,
                         a = pr$a,
                         b = pr$b,
                         c = pr$c,
                         prior_var_phi = pr$prior_var_phi,
                         prior_delta_c = pr$prior_delta_c,
                         prior_delta_d = pr$prior_delta_d,
                         print_progress = print_progress,
                         user_seed = user_seed)

    result <- DetectCpObj(data = data,
                          n_iterations = n_iterations,
                          n_burnin = n_burnin,
                          orders = out$orders,
                          time = out$time,
                          entropy_MCMC = coda::mcmc(out$entropy),
                          lkl_MCMC = coda::mcmc(out$lkl),
                          phi_MCMC = coda::mcmc(out$phi_MCMC),
                          sigma_MCMC = coda::mcmc(out$sigma_MCMC),
                          delta_MCMC = coda::mcmc(out$delta_MCMC),
                          kernel_ts = TRUE,
                          kernel_epi = FALSE,
                          univariate_ts = TRUE)

  } else if(kernel == "ts"){

    rp <- resolve_params(params, "multi", d = nrow(data))
    report_sampler_settings(rp$sampler, "detect_cp")
    pr <- rp$values

    out <- detect_cp_multi(data = as.matrix(data),
                           n_iterations = n_iterations,
                           q = q,
                           m_0 = pr$m_0,
                           k_0 = pr$k_0,
                           nu_0 = pr$nu_0,
                           S_0 = pr$S_0,
                           prior_delta_c = pr$prior_delta_c,
                           prior_delta_d = pr$prior_delta_d,
                           prior_var_phi = pr$prior_var_phi,
                           print_progress = print_progress,
                           user_seed = user_seed)

    result <- DetectCpObj(data = data,
                          n_iterations = n_iterations,
                          n_burnin = n_burnin,
                          orders = out$orders,
                          time = out$time,
                          entropy_MCMC = coda::mcmc(out$entropy),
                          lkl_MCMC = coda::mcmc(out$lkl),
                          phi_MCMC = coda::mcmc(out$phi_MCMC),
                          sigma_MCMC = coda::mcmc(out$sigma_MCMC),
                          delta_MCMC = coda::mcmc(out$delta_MCMC),
                          kernel_ts = TRUE,
                          kernel_epi = FALSE,
                          univariate_ts = FALSE)

  } else {

    rp <- resolve_params(params, "epi")
    report_sampler_settings(rp$sampler, "detect_cp")
    pr <- rp$values

    out <- detect_cp_epi(data = data,
                         n_iterations = n_iterations,
                         q = q,
                         M = pr$M,
                         xi = pr$xi,
                         a0 = pr$a0,
                         b0 = pr$b0,
                         I0_var = pr$I0_var,
                         print_progress = print_progress,
                         user_seed = user_seed)

    result <- DetectCpObj(data = data,
                          n_iterations = n_iterations,
                          n_burnin = n_burnin,
                          orders = out$orders,
                          time = out$time,
                          entropy_MCMC = coda::mcmc(out$entropy),
                          lkl_MCMC = coda::mcmc(out$lkl),
                          I0_MCMC = coda::mcmc(out$I0_MCMC),
                          kernel_ts = FALSE,
                          kernel_epi = TRUE,
                          univariate_ts = FALSE)

  }

  return(result)

}
