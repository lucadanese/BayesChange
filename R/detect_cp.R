#' @name detect_cp
#' @export detect_cp
#'
#' @title Detect change points on time series.
#' @description The \code{detect_cp} function detect change points on univariate and multivariate time series.
#'
#'
#' @param data a vector or a matrix. If a vector the algorithm for
#' univariate time series is used. If a matrix, where rows are the observations
#' and columns are the times, then the algorithm for multivariate time series is used.
#'
#' @param n_iterations number of MCMC iterations.
#' @param n_burnin number of iterations that must be excluded when computing the posterior estimate.
#' @param print_progress If TRUE (default) print the progress bar.
#' @param user_seed seed for random distribution generation.
#'
#'
#' @param params a list of parameters:
#'
#' If the time series is univariate the following must be specified:
#'
#' \itemize{
#'   \item \code{q} probability of performing a split at each iteration.
#'   \item \code{phi} parameter \eqn{\phi} of the integrated likelihood function.
#'   \item \code{a}, \code{b}, \code{c} parameters of the Normal-Gamma prior for \eqn{\mu} and \eqn{\lambda}.
#'   \item \code{par_theta_c}, \code{par_theta_d} parameters of the shifted Gamma prior for \eqn{\theta}.
#' }
#'
#' If the time series is multivariate the following must be specified:
#'
#' \itemize{
#'   \item \code{q} probability of performing a split at each iteration.
#'   \item \code{k_0}, \code{nu_0}, \code{phi_0}, \code{m_0} parameters for the Normal-Inverse-Wishart prior for \eqn{(\mu,\lambda)}.
#'   \item \code{par_theta_c}, \code{par_theta_d} parameters for the shifted Gamma prior for \eqn{\theta}.
#'   \item \code{prior_var_gamma} parameters for the Gamma prior for \eqn{\gamma}.
#'   \item \code{print_progress} If TRUE (default) print the progress bar.
#'   \item \code{user_seed} seed for random distribution generation.
#' }
#'
#' @return A \code{DetectCpObj} class object containing
#'
#' \itemize{
#'   \item \code{$data} vector or matrix with the data.
#'   \item \code{$n_iterations} number of iterations.
#'   \item \code{$n_burnin} number of burn-in iterations.
#'   \item \code{$orders} matrix where each entries is the assignment of the realization to a block. Rows are the iterations and columns the times.
#'   \item \code{$time} computational time.
#'   \item \code{$gammaMCMC} traceplot for \eqn{\gamma}.
#'   \item \code{$gamma_MCMC_01} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\gamma} was accepted, \eqn{0} otherwise.
#'   \item \code{$sigma_MCMC} traceplot for \eqn{\sigma}.
#'   \item \code{$sigma_MCMC_01} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\sigma} was accepted, \eqn{0} otherwise.
#'   \item \code{$theta_MCMC} traceplot for \eqn{\theta}.
#'   \item \code{$univariate_ts} TRUE if data is an univariate time series, FALSE if it is a multivariate time series.
#' }
#'
#' @examples
#'
#' ## Univariate time series
#'
#' data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))
#'
#'
#' out <- detect_cp(data = data_vec, n_iterations = 2500, n_burnin = 500,
#'                  params = list(q = 0.25, phi = 0.1, a = 1, b = 1, c = 0.1))
#'
#' print(out)
#'
#' ## Multivariate time series
#'
#' data_mat <- matrix(NA, nrow = 3, ncol = 100)
#'
#' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
#' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#'
#'
#' out <- detect_cp(data = data_mat, n_iterations = 2500, n_burnin = 500,
#'                  params = list(q = 0.25, k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
#'                                par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1))
#'
#' print(out)
#'
#' @references
#'
#' Martínez, A. F., & Mena, R. H. (2014). On a Nonparametric Change Point Detection Model in Markovian Regimes. \emph{Bayesian Analysis}, 9(4), 823–858. \doi{10.1214/14-BA878}
#'
#' Corradin, R., Danese, L., & Ongaro, A. (2022). Bayesian nonparametric change point detection for multivariate time series with missing observations. \emph{International Journal of Approximate Reasoning}, 143, 26--43. \doi{10.1016/j.ijar.2021.12.019}
#'
#'
detect_cp <- function(data,
                      n_iterations,
                      n_burnin = 0,
                      params = list(),
                      print_progress = TRUE,
                      user_seed = 1234){

  if((!is.vector(data)) && (!is.matrix(data))) stop("data must be a vector or a matrix")
  if(is.null(n_iterations)) stop("missing number of iterations")
  if(n_iterations < 1) stop("number of iterations must be positive and > 1")
  if(n_burnin < 0) stop("number of burn in iterations must be positive and > 1")


  if(is.vector(data)){

    if((!is.null(params$q)) && ((params$q >= 1) | (params$q <= 0))) stop("params$q must be in (0,1)")
    if((!is.null(params$phi)) && ((params$phi >= 1) | (params$phi <= 0))) stop("params$phi must be in (0,1)")
    if((!is.null(params$a)) && (params$a < 0)) stop("params$a must be positive")
    if((!is.null(params$b)) && (params$b < 0)) stop("params$b must be positive")
    if(!(is.null(params$c)) && (params$c < 0)) stop("params$c must be positive")
    if(!is.null(params$par_theta_c) && (params$par_theta_c < 0)) stop("params$par_theta_c must be positive")
    if(!is.null(params$par_theta_d) && (params$par_theta_d < 0)) stop("params$par_theta_d must be positive")
    if((!is.null(params$params) && !is.list(params))) stop("params must be a list")
    #if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
    #if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

    # substitute missing parameters with default

    q_input = ifelse(is.null(params$q), 0.5, params$q)
    n_burnin_input = n_burnin
    phi_input = ifelse(is.null(params$phi), 0.1, params$phi)
    a_input = ifelse(is.null(params$a), 1, params$a)
    b_input = ifelse(is.null(params$b), 1, params$b)
    c_input = ifelse(is.null(params$c), 0.1, params$c)
    par_theta_c_input = ifelse(is.null(params$par_theta_c), 1, params$par_theta_c)
    par_theta_d_input = ifelse(is.null(params$par_theta_d), 1, params$par_theta_d)
    #print_progress_input = ifelse(is.null(print_progress), TRUE, print_progress)
    #user_seed_input = ifelse(is.null(user_seed), 1234, user_seed)
    print_progress_input = print_progress
    user_seed_input = user_seed


    #

    data_input = data
    n_iterations_input = n_iterations

    out <- detect_cp_uni(data = data_input,
                         n_iterations = n_iterations_input,
                         q = q_input,
                         phi = phi_input,
                         a = a_input,
                         b = b_input,
                         c = c_input,
                         par_theta_c = par_theta_c_input,
                         par_theta_d = par_theta_d_input,
                         print_progress = print_progress_input,
                         user_seed = user_seed_input)

    # save output

    result <- DetectCpObj(data = data_input,
                          n_iterations = n_iterations_input,
                          n_burnin = n_burnin_input,
                          orders = out$orders,
                          time = out$time,
                          sigma_MCMC = out$sigma_MCMC,
                          sigma_MCMC_01 = out$sigma_MCMC_01,
                          theta_MCMC = out$theta_MCMC,
                          univariate_ts = TRUE)

  }

  if(is.matrix(data)){

    if((!is.null(params$q)) && ((params$q >= 1) | (params$q <= 0))) stop("params$q must be in (0,1)")
    if((!is.null(params$k_0)) && (params$k_0 < 0)) stop("params$k_0 must be positive")
    if((!is.null(params$nu_0)) && (params$nu_0 < 0)) stop("params$nu_0 must be positive")
    if((!is.null(params$phi_0)) && nrow(params$phi_0) != ncol(params$phi_0)) stop("number of rows and columns must be the same in params$phi_0")
    if((!is.null(params$m_0)) && (length(params$m_0) != nrow(data))) stop("number of elements in params$m_0 must equal the number of observations")
    if((!is.null(params$par_theta_c)) && (params$par_theta_c < 0)) stop("params$par_theta_c must be positive")
    if((!is.null(params$par_theta_d)) && (params$par_theta_d < 0)) stop("params$par_theta_d must be positive")
    if((!is.null(params$prior_var_gamma)) && (params$prior_var_gamma < 0)) stop("params$prior_var_gamma must be positive")
    if((!is.null(params$params) && !is.list(params))) stop("params must be a list")
    if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
    if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

    # substitute missing parameters with default

    q_input = ifelse(is.null(params$q), 0.5, params$q)
    n_burnin_input = n_burnin
    k_0_input = ifelse(is.null(params$k_0), 0.5, params$k_0)
    nu_0_input = ifelse(is.null(params$nu_0), nrow(data)+1, params$nu_0)
    par_theta_c_input = ifelse(is.null(params$par_theta_c), 1, params$par_theta_c)
    par_theta_d_input = ifelse(is.null(params$par_theta_d), 1, params$par_theta_d)
    prior_var_gamma_input = ifelse(is.null(params$prior_var_gamma), 0.1, params$prior_var_gamma)
    #print_progress_input = ifelse(is.null(print_progress), TRUE, print_progress)
    #user_seed_input = ifelse(is.null(user_seed), 1234, user_seed)
    print_progress_input = print_progress
    user_seed_input = user_seed

    # with object matrix ifelse does not work
    if(is.null(params$phi_0)){phi_0_input = diag(0.1, nrow(data), nrow(data))} else{phi_0_input = params$phi_0}
    if(is.null(params$m_0)){m_0_input = rep(0, nrow(data))} else{m_0_input = params$m_0}
    #

    #

    data_input = data
    n_iterations_input = n_iterations
    n_burnin_input = n_burnin

    out <- detect_cp_multi(data = as.matrix(data_input),
                           n_iterations = n_iterations_input,
                           q = q_input,
                           k_0 = k_0_input,
                           nu_0 = nu_0_input,
                           phi_0 = phi_0_input,
                           m_0 = m_0_input,
                           par_theta_c = par_theta_c_input,
                           par_theta_d = par_theta_d_input,
                           prior_var_gamma = prior_var_gamma_input,
                           print_progress = print_progress_input,
                           user_seed = user_seed_input)

    # save output

    result <- DetectCpObj(data = data_input,
                          n_iterations = n_iterations_input,
                          n_burnin = n_burnin_input,
                          orders = out$orders,
                          time = out$time,
                          gamma_MCMC = out$gamma_MCMC,
                          gamma_MCMC_01 = out$gamma_MCMC_01,
                          sigma_MCMC = out$sigma_MCMC,
                          sigma_MCMC_01 = out$sigma_MCMC_01,
                          theta_MCMC = out$theta_MCMC,
                          univariate_ts = FALSE)

  }

  return(result)


}
