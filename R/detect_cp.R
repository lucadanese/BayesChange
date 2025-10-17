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
#' @param user_seed seed for random distribution generation.
#'
#'
#' @param params a list of parameters:
#'
#' If data is an univariate time series the following must be specified:
#'
#' \itemize{
#'   \item \code{a}, \code{b}, \code{c} parameters of the Normal-Gamma prior for \eqn{\mu} and \eqn{\lambda}.
#'   \item \code{prior_var_phi} variance for the proposal in the \eqn{N(0,\sigma^2_\phi)} posterior estimate of \eqn{\delta}.
#'   \item \code{prior_delta_c} parameter of the shifted Gamma prior of \eqn{\delta}.
#'   \item \code{prior_delta_d} parameter of the shifted Gamma prior of \eqn{\delta}.
#' }
#'
#' If the time series is multivariate the following must be specified:
#'
#' \itemize{
#'   \item \code{m_0}, \code{k_0}, \code{nu_0}, \code{S_0} parameters for the Normal-Inverse-Wishart prior for \eqn{(\mu,\lambda)}.
#'   \item \code{prior_var_phi} variance for the proposal in the \eqn{N(0,\sigma^2_\phi)} posterior estimate of \eqn{\delta}.
#'   \item \code{prior_delta_c} parameter of the shifted Gamma prior of \eqn{\delta}.
#'   \item \code{prior_delta_d} parameter of the shifted Gamma prior of \eqn{\delta}.
#' }
#'
#' If data are epidemic diffusions:
#'
#' \itemize{
#'   \item \code{M} number of Monte Carlo iterations when computing the likelihood of the epidemic diffusion.
#'   \item \code{xi} recovery rate fixed constant for each population at each time.
#'   \item \code{a0},\code{b0} parameters for the computation of the integrated likelihood of the epidemic diffusions.
#'   \item \code{I0_var} variance for the Metropolis-Hastings estimation of the proportion of infected at time 0.
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
#'   \item \code{$phi_MCMC} traceplot for \eqn{\gamma}.
#'   \item \code{$phi_MCMC_01} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\gamma} was accepted, \eqn{0} otherwise.
#'   \item \code{$sigma_MCMC} traceplot for \eqn{\sigma}.
#'   \item \code{$sigma_MCMC_01} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\sigma} was accepted, \eqn{0} otherwise.
#'   \item \code{$delta_MCMC} traceplot for \eqn{\delta}.
#'   \item \code{I0_MCMC} traceplot for \eqn{I_0}.
#'   \item \code{I0_MCMC_01} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{I_0} was accepted, \eqn{0} otherwise.
#'   \item \code{kernel_ts} if TRUE data are time series.
#'   \item \code{kernel_epi} if TRUE data are epidemic diffusions.
#'   \item \code{$univariate_ts} TRUE if data is an univariate time series, FALSE if it is a multivariate time series.
#' }
#'
#' @examples
#'
#' ## Univariate time series
#'
#' data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))
#'
#' out <- detect_cp(data = data_vec, n_iterations = 2500, n_burnin = 500,
#'                  params = list(a = 1, b = 1, c = 0.1), kernel = "ts")
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
#' out <- detect_cp(data = data_mat, n_iterations = 2500, n_burnin = 500,
#'                  params = list(m_0 = rep(0,3), k_0 = 0.25, nu_0 = 4,
#'                                S_0 = diag(1,3,3)), kernel = "ts")
#'
#' print(out)
#'
#' \donttest{
#' ## Epidemic diffusions
#'
#' data_mat <- matrix(NA, nrow = 100, ncol = 1)
#'
#' betas <- c(rep(0.45, 25),rep(0.14,75))
#'
#' inf_times <- sim_epi_data(10000, 10, 100, betas, 1/8)
#'
#' inf_times_vec <- rep(0,100)
#' names(inf_times_vec) <- as.character(1:100)
#'
#' for(j in 1:100){
#'  if(as.character(j) %in% names(table(floor(inf_times)))){
#'    inf_times_vec[j] = table(floor(inf_times))[which(names(table(floor(inf_times))) == j)]
#'  }
#' }
#'
#' data_mat[,1] <- inf_times_vec
#'
#' out <- detect_cp(data = data_mat, n_iterations = 500, n_burnin = 100,
#'                  params = list(M = 250, xi = 1/8, a0 = 40, b0 = 10), kernel = "epi")
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
#' Corradin, R., Danese, L., KhudaBukhsh, W. R., & Ongaro, A. (2024). \emph{Model-based clustering of time-dependent observations with common structural changes}. arXiv preprint arXiv:2410.09552.
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

  if(kernel == "ts"){

    if(is.vector(data)){

      if((!is.null(params$a)) && (params$a < 0)) stop("params$a must be positive")
      if((!is.null(params$b)) && (params$b < 0)) stop("params$b must be positive")
      if(!(is.null(params$c)) && (params$c < 0)) stop("params$c must be positive")
      if((!is.null(params$params) && !is.list(params))) stop("params must be a list")
      if((!is.null(params$prior_var_phi)) && (params$prior_var_phi <= 0)) stop("prior_var_phi must be positive")
      if(!is.null(params$prior_delta_c) && (params$prior_delta_c < 0)) stop("prior_delta_c must be positive")
      if(!is.null(params$prior_delta_d) && (params$prior_delta_d < 0)) stop("prior_delta_d must be positive")

      data_input = data
      n_iterations_input = n_iterations
      n_burnin_input = n_burnin
      q_input = q
      print_progress_input = print_progress
      user_seed_input = user_seed

      # substitute missing parameters with default
      a_input = ifelse(is.null(params$a), 1, params$a)
      b_input = ifelse(is.null(params$b), 1, params$b)
      c_input = ifelse(is.null(params$c), 0.1, params$c)
      prior_var_phi_input = ifelse(is.null(params$prior_var_phi), 0.1, params$prior_var_phi)
      prior_delta_c_input = ifelse(is.null(params$prior_delta_c), 1, params$prior_delta_c)
      prior_delta_d_input = ifelse(is.null(params$prior_delta_d), 1, params$prior_delta_d)
      #

      out <- detect_cp_uni(data = data_input,
                           n_iterations = n_iterations_input,
                           q = q_input,
                           a = a_input,
                           b = b_input,
                           c = c_input,
                           prior_var_phi = prior_var_phi_input,
                           prior_delta_c = prior_delta_c_input,
                           prior_delta_d = prior_delta_d_input,
                           print_progress = print_progress_input,
                           user_seed = user_seed_input)

      # save output

      result <- DetectCpObj(data = data_input,
                            n_iterations = n_iterations_input,
                            n_burnin = n_burnin_input,
                            orders = out$orders,
                            time = out$time,
                            phi_MCMC = out$phi_MCMC,
                            phi_MCMC_01 = out$phi_MCMC_01,
                            sigma_MCMC = out$sigma_MCMC,
                            sigma_MCMC_01 = out$sigma_MCMC_01,
                            delta_MCMC = out$delta_MCMC,
                            kernel_ts = TRUE,
                            kernel_epi = FALSE,
                            univariate_ts = TRUE)

    }

    if(is.matrix(data)){

      if((!is.null(params$k_0)) && (params$k_0 < 0)) stop("params$k_0 must be positive")
      if((!is.null(params$nu_0)) && (params$nu_0 < 0)) stop("params$nu_0 must be positive")
      if((!is.null(params$S_0)) && nrow(params$S_0) != ncol(params$S_0)) stop("number of rows and columns must be the same in params$S_0")
      if((!is.null(params$m_0)) && (length(params$m_0) != nrow(data))) stop("number of elements in params$m_0 must equal the number of observations")

      # substitute missing parameters with default
      if(is.null(params$m_0)){m_0_input = rep(0, nrow(data))} else{m_0_input = params$m_0}
      k_0_input = ifelse(is.null(params$k_0), 0.5, params$k_0)
      nu_0_input = ifelse(is.null(params$nu_0), nrow(data)+1, params$nu_0)
      if(is.null(params$S_0)){S_0_input = diag(0.1, nrow(data), nrow(data))} else{S_0_input = params$S_0}
      prior_var_phi_input = ifelse(is.null(params$prior_var_phi), 0.1, params$prior_var_phi)
      prior_delta_c_input = ifelse(is.null(params$prior_delta_c), 1, params$prior_delta_c)
      prior_delta_d_input = ifelse(is.null(params$prior_delta_d), 1, params$prior_delta_d)
      #

      data_input = data
      n_iterations_input = n_iterations
      n_burnin_input = n_burnin
      q_input = q
      print_progress_input = print_progress
      user_seed_input = user_seed

      out <- detect_cp_multi(data = as.matrix(data_input),
                             n_iterations = n_iterations_input,
                             q = q_input,
                             m_0 = m_0_input,
                             k_0 = k_0_input,
                             nu_0 = nu_0_input,
                             S_0 = S_0_input,
                             prior_delta_c = prior_delta_c_input,
                             prior_delta_d = prior_delta_d_input,
                             prior_var_phi = prior_var_phi_input,
                             print_progress = print_progress_input,
                             user_seed = user_seed_input)

      # save output

      result <- DetectCpObj(data = data_input,
                            n_iterations = n_iterations_input,
                            n_burnin = n_burnin_input,
                            orders = out$orders,
                            time = out$time,
                            phi_MCMC = out$phi_MCMC,
                            phi_MCMC_01 = out$phi_MCMC_01,
                            sigma_MCMC = out$sigma_MCMC,
                            sigma_MCMC_01 = out$sigma_MCMC_01,
                            delta_MCMC = out$delta_MCMC,
                            kernel_ts = TRUE,
                            kernel_epi = FALSE,
                            univariate_ts = FALSE)

    }

  }

  if(kernel == "epi"){

    if((!is.null(params$M)) && (params$M < 1)) stop("params$M must be at least equal to 1")
    if((!is.null(params$xi)) && ((params$xi <= 0) | (params$xi >= 1))) stop("params$xi must be in (0,1)")
    if((!is.null(params$alpha_SM)) && (params$alpha_SM <= 0)) stop("params$alpha_SM must be positive")
    if((!is.null(params$a0)) && (params$a0 <= 0)) stop("params$a0 must be positive")
    if((!is.null(params$b0)) && (params$b0 <= 0)) stop("params$b0 must be positive")
    if((!is.null(params$I0_var)) && ((params$I0_var <= 0) | (params$I0_var >= 1))) stop("params$I0_var must be in (0,1)")
    if((!is.null(params$avg_blk) && (params$avg_blk < 0))) stop("params$avg_blk must be positive")
    if((!is.null(params) && !is.list(params))) stop("params must be a list")

    data_input = data
    n_iterations_input = n_iterations
    n_burnin_input = n_burnin
    q_input = q
    print_progress_input = print_progress
    user_seed_input = user_seed

    # substitute missing parameters with default
    M_input = ifelse(is.null(params$M), 250, params$M)
    xi_input = ifelse(is.null(params$xi), 1/8, params$xi)
    a0_input = ifelse(is.null(params$a0), 4, params$a0)
    b0_input = ifelse(is.null(params$b0), 10, params$b0)
    I0_var_input = ifelse(is.null(params$I0_var), 0.01, params$I0_var)
    #

    out <- detect_cp_epi(data = data_input,
                         n_iterations = n_iterations_input,
                         q = q_input,
                         M = M_input,
                         xi = xi_input,
                         a0 = a0_input,
                         b0 = b0_input,
                         I0_var = I0_var_input,
                         print_progress = print_progress_input,
                         user_seed = user_seed_input)

    result <- DetectCpObj(data = data_input,
                          n_iterations = n_iterations_input,
                          n_burnin = n_burnin_input,
                          orders = out$orders,
                          time = out$time,
                          I0_MCMC = out$I0_MCMC,
                          I0_MCMC_01 = out$I0_MCMC_01,
                          kernel_ts = FALSE,
                          kernel_epi = TRUE,
                          univariate_ts = FALSE)

  }



  return(result)


}
