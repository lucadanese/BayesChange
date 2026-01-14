#' @name clust_cp
#' @export clust_cp
#'
#' @title Clustering time dependent observations with common change points.
#' @description The \code{clust_cp} function cluster observations with common change points. Data can be time series or epidemic diffusions.
#'
#' @param data a matrix or an array If a matrix the algorithm for
#' univariate time series is used, where each row is a time series. If an array, the algorithm is run for multivariate time series. Each slice of the array is a matrix where the rows are the dimensions of the time series.
#' @param n_iterations number of MCMC iterations.
#' @param n_burnin number of iterations that must be excluded when computing the posterior estimate.
#' @param alpha_SM \eqn{\alpha} for the split-merge main algorithm.
#' @param B number of orders for the normalization constant.
#' @param L number of split-merge steps for the proposal step.
#' @param q probability of a split in the split-merge proposal and acceleration step.
#' @param print_progress If TRUE (default) print the progress bar.
#' @param user_seed seed for random distribution generation.
#' @param kernel can be "ts" if data are time series or "epi" if data are epidemic diffusions.
#'
#' @param params a list of parameters:
#'
#' If the time series is univariate the following must be specified:
#'
#' \itemize{
#'   \item \code{a},\code{b},\code{c} parameters of the integrated likelihood.
#'   \item \code{phi} correlation parameter in the likelihood.
#' }
#'
#' If the time series is multivariate the following must be specified:
#'
#' \itemize{
#'   \item \code{k_0}, \code{nu_0}, \code{S_0}, \code{m_0} parameters of the integrated likelihood.
#'   \item \code{phi} correlation parameter in the likelihood.
#' }
#'
#' If data are epidemic diffusions:
#'
#' \itemize{
#'   \item \code{M} number of Monte Carlo iterations when computing the likelihood of the epidemic diffusion.
#'   \item \code{xi} recovery rate fixed constant for each population at each time.
#'   \item \code{a0}, \code{b0} parameters for the computation of the integrated likelihood of the epidemic diffusions.
#'   \item \code{I0_var} variance for the Metropolis-Hastings estimation of the proportion of infected at time 0.
#'   \item \code{avg_blk} prior average number of change points for each order.
#' }
#'
#' @return A \code{ClustCpObj} class object containing
#'
#' \itemize{
#'   \item \code{$data} Vector or matrix containing the data.
#'   \item \code{$n_iterations} Total number of MCMC iterations.
#'   \item \code{$n_burnin} Number of burn-in iterations.
#'   \item \code{$clust} A matrix where each row corresponds to the cluster assignment from each iteration.
#'   \item \code{$orders} A multidimensional array where each slice is a matrix representing the latent order at each iteration.
#'   \item \code{$time} Total computational time (in seconds).
#'   \item \code{$entropy_MCMC} A \code{coda::mcmc} object containing the MCMC samples of the entropy.
#'   \item \code{$lkl} A \code{coda::mcmc} object containing the log-likelihood evaluated at each iteration.
#'   \item \code{$norm_vec} A vector containing the normalization constants computed at the beginning of the algorithm.
#'   \item \code{$I0_MCMC} A \code{coda::mcmc} object containing the MCMC trace of the initial infection proportion \eqn{I_0}.
#'   \item \code{$kernel_ts} TRUE if the kernel used corresponds to time series data.
#'   \item \code{$kernel_epi} TRUE if the kernel used corresponds to epidemic diffusion data.
#'   \item \code{$univariate_ts} TRUE if the data represent a univariate time series, FALSE if multivariate.
#' }
#'
#'
#' @examples
#'
#'\donttest{
#' ## Univariate time series
#'
#' data("stock_uni")
#'
#' params_uni <- list(a = 1,
#'                    b = 1,
#'                    c = 1,
#'                    phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:5,], n_iterations = 2000, n_burnin = 500,
#'                 L = 1, q = 0.5, B = 1000, params = params_uni, kernel = "ts")
#'
#' print(out)
#'
#' ## Multivariate time series
#'
#' data("stock_multi")
#'
#' params_multi <- list(m_0 = rep(0,2),
#'                      k_0 = 1,
#'                      nu_0 = 10,
#'                      S_0 = diag(1,2,2),
#'                      phi = 0.1)
#'
#' out <- clust_cp(data = stock_multi[,,1:5], n_iterations = 2000, n_burnin = 500,
#'                 L = 1, B = 1000, params = params_multi, kernel = "ts")
#'
#' print(out)
#'
#' ## Epidemic diffusions
#'
#' data("epi_synthetic_multi")
#'
#' params_epi <- list(M = 100, xi = 1/8,
#'                    alpha_SM = 1,
#'                    a0 = 4,
#'                    b0 = 10,
#'                    I0_var = 0.1,
#'                    avg_blk = 2)
#'
#' out <- clust_cp(epi_synthetic_multi, n_iterations = 2000, n_burnin = 500,
#'                 L = 1, B = 1000, params = params_epi, kernel = "epi")
#'
#' print(out)
#'
#' }
#'
#' @references
#'
#' Corradin, R., Danese, L., KhudaBukhsh, W. R., & Ongaro, A. (2026). Model-based clustering of time-dependent observations with common structural changes. \emph{Statistics and Computing}. \doi{10.1007/s11222-025-10756-x}
#'
clust_cp <- function(data,
                     n_iterations,
                     n_burnin = 0,
                     params = list(),
                     alpha_SM = 1,
                     B = 1000,
                     L = 1,
                     q = 0.5,
                     kernel,
                     print_progress = TRUE,
                     user_seed = 1234){

  if((!is.matrix(data)) && (!is.array(data))) stop("data must be a matrix or an array")
  if(n_iterations < 1) stop("number of iterations must be positive and > 1")
  if(n_burnin < 0) stop("number of burn in iterations must be positive and > 1")
  if(alpha_SM < 0) stop("alpha_SM must be positive")
  if((!is.null(B)) && (B < 1)) stop("B must be at least equal to 1")
  if((!is.null(L)) && (L < 1)) stop("L must be at least equal to 1")
  if((!is.null(q)) && ((q >= 1) | (q <= 0))) stop("params$q must be in (0,1)")
  if((is.null(kernel) | !(kernel %in% c("ts", "epi")))) stop("kernel must be 'ts' or 'epi'")
  if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
  if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

  if(kernel == "ts"){

    if(is.matrix(data)){

      if((!is.null(params$phi)) && ((params$phi > 1) | (params$phi < 0))) stop("params$phi must be in (0,1)")
      if((!is.null(params$a)) && (params$a <= 0)) stop("params$a must be positive")
      if((!is.null(params$b)) && (params$b <= 0)) stop("params$b must be positive")
      if((!is.null(params$c)) && (params$c <= 0)) stop("params$c must be positive")
      if(!is.null(params) && !is.list(params)) stop("params must be a list")


      data_input = data
      n_iterations_input = n_iterations
      alpha_SM_input = alpha_SM
      n_burnin_input = n_burnin
      B_input = B
      L_input = L
      print_progress_input = print_progress
      user_seed_input = user_seed

      # substitute missing parameters with default
      phi_input = ifelse(is.null(params$phi), 0.1, params$phi)
      a_input = ifelse(is.null(params$a), 1, params$a)
      b_input = ifelse(is.null(params$b), 1, params$b)
      c_input = ifelse(is.null(params$c), 0.1, params$c)
      q_input = ifelse(is.null(params$q), 0.5, params$q)
      #

      out <- clust_cp_uni(data = data_input,
                          n_iterations = n_iterations_input,
                          B = B_input, L = L_input,
                          phi = phi_input, a = a_input, b = b_input,
                          c = c_input, q = q_input,
                          alpha_SM = alpha_SM_input,
                          print_progress = print_progress_input, user_seed = user_seed_input)

      # create coda objects
      entropy_MCMC <- coda::mcmc(out$entropy)
      lkl_MCMC <- coda::mcmc(out$lkl)

      result <- ClustCpObj(data = data_input,
                           n_iterations = n_iterations_input,
                           n_burnin = n_burnin_input,
                           clust = out$clust,
                           orders = out$orders,
                           time = out$time,
                           norm_vec = out$norm_vec,
                           entropy_MCMC = entropy_MCMC,
                           lkl_MCMC = lkl_MCMC,
                           kernel_ts = TRUE,
                           kernel_epi = FALSE,
                           univariate_ts = TRUE)

    } else if(is.array(data)){

      if((!is.null(params$k_0)) && (params$k_0 < 0)) stop("params$k_0 must be positive")
      if((!is.null(params$nu_0)) && (params$nu_0 < 0)) stop("params$nu_0 must be positive")
      if((!is.null(params$S_0)) && nrow(params$S_0) != ncol(params$S_0)) stop("number of rows and columns must be the same in params$S_0")
      if((!is.null(params$m_0)) && (length(params$m_0) != nrow(data))) stop("number of elements in params$m_0 must equal the number of observations")
      if((!is.null(params$q)) && ((params$q >= 1) | (params$q <= 0))) stop("params$q must be in (0,1)")
      if((!is.null(params$params) && !is.list(params))) stop("params must be a list")

      data_input = data
      n_iterations_input = n_iterations
      n_burnin_input = n_burnin
      alpha_SM_input = alpha_SM
      B_input = B
      L_input = L
      q_input = q
      print_progress_input = print_progress
      user_seed_input = user_seed

      # substitute missing parameters with default
      if(is.null(params$m_0)){m_0_input = rep(0, nrow(data))} else{m_0_input = params$m_0}
      k_0_input = ifelse(is.null(params$k_0), 0.5, params$k_0)
      phi_input = ifelse(is.null(params$phi), 0.1, params$phi)
      if(is.null(params$S_0)){S_0_input = diag(0.1, nrow(data), nrow(data))} else{S_0_input = params$S_0}
      nu_0_input = ifelse(is.null(params$nu_0), nrow(data)+1, params$nu_0)
      #

      out <- clust_cp_multi(data = data_input, n_iterations = n_iterations_input,
                            B = B_input, L = L_input, phi = phi_input,
                            k_0 = k_0_input, nu_0 = nu_0_input,
                            S_0 = S_0_input, m_0 = m_0_input, q = q_input,
                            alpha_SM = alpha_SM_input,
                            print_progress = print_progress_input,
                            user_seed = user_seed_input)

      # create coda objects
      entropy_MCMC <- coda::mcmc(out$entropy)
      lkl_MCMC <- coda::mcmc(out$lkl)

      result <- ClustCpObj(data = data_input,
                           n_iterations = n_iterations_input,
                           n_burnin = n_burnin_input,
                           clust = out$clust,
                           orders = out$orders,
                           time = out$time,
                           norm_vec = out$norm_vec,
                           entropy_MCMC = entropy_MCMC,
                           lkl_MCMC = lkl_MCMC,
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
    alpha_SM_input = alpha_SM
    B_input = B
    L_input = L
    q_input = q
    print_progress_input = print_progress
    user_seed_input = user_seed

    # substitute missing parameters with default
    M_input = ifelse(is.null(params$M), 250, params$M)
    xi_input = ifelse(is.null(params$xi), 1/8, params$xi)
    a0_input = ifelse(is.null(params$a0), 4, params$a0)
    b0_input = ifelse(is.null(params$b0), 10, params$b0)
    I0_var_input = ifelse(is.null(params$I0_var), 0.01, params$I0_var)
    avg_blk_input = ifelse(is.null(params$avg_blk), 2, params$avg_blk)
    #

    out <- clust_cp_epi(data = data_input, n_iterations = n_iterations_input,
                        M = M_input, B = B_input, L = L_input,
                        xi = xi_input,alpha_SM = alpha_SM_input, q = q_input,
                        a0 = a0_input, b0 = b0_input, I0_var = I0_var_input,
                        avg_blk = avg_blk_input, print_progress = print_progress_input,
                        user_seed = user_seed_input)

    # create coda objects
    entropy_MCMC <- coda::mcmc(out$entropy)
    I0_MCMC <- coda::mcmc(out$I0_MCMC)
    lkl_MCMC <- coda::mcmc(out$llik)

    result <- ClustCpObj(data = data_input,
                         n_iterations = n_iterations_input,
                         n_burnin = n_burnin_input,
                         clust = out$clust,
                         orders = out$orders,
                         time = out$time,
                         norm_vec = out$norm_vec,
                         entropy_MCMC = entropy_MCMC,
                         lkl_MCMC = lkl_MCMC,
                         I0_MCMC = I0_MCMC,
                         kernel_ts = FALSE,
                         kernel_epi = TRUE,
                         univariate_ts = TRUE)

  }

  return(result)

}
