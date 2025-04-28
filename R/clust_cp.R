#' @name clust_cp
#' @export clust_cp
#'
#' @title Clustering time dependent observations with common change points.
#' @description The \code{clust_cp} function cluster observations with common change points. Data can be time series or survival functions.
#'
#' @param data a matrix or an array If a matrix the algorithm for
#' univariate time series is used, where each row is a time series. If an array, the algorithm is run for multivariate time series. Each slice of the array is a matrix where the rows are the dimensions of the time series.
#' @param n_iterations number of MCMC iterations.
#' @param n_burnin number of iterations that must be excluded when computing the posterior estimate.
#' @param B number of orders for the normalization constant.
#' @param L number of split-merge steps for the proposal step.
#' @param q probability of a split in the split-merge proposal and acceleration step.
#' @param print_progress If TRUE (default) print the progress bar.
#' @param alpha_SM \eqn{\alpha} for the split-merge main algorithm.
#' @param user_seed seed for random distribution generation.
#' @param kernel can be "ts" if data are time series or "epi" if data are survival functions.
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
#'   \item \code{k_0},\code{nu_0},\code{S_0},\code{m_0} parameters of the integrated likelihood.
#'   \item \code{phi} correlation parameter in the likelihood.
#' }
#'
#' If data are survival functions:
#'
#' \itemize{
#'   \item \code{M} number of Monte Carlo iterations when computing the likelihood of the survival function.
#'   \item \code{xi} recovery rate fixed constant for each population at each time.
#'   \item \code{a0},\code{b0} parameters for the computation of the integrated likelihood of the survival functions.
#'   \item \code{I0_var} variance for the Metropolis-Hastings estimation of the proportion of infected at time 0.
#'   \item \code{p} prior average number of change points for each order.
#' }
#'
#' @return A \code{ClustCpObj} class object containing
#'
#' \itemize{
#'   \item \code{$data} vector or matrix with the data.
#'   \item \code{$n_iterations} number of iterations.
#'   \item \code{$n_burnin} number of burn-in iterations.
#'   \item \code{$clust} a matrix where each row corresponds to the output cluster of the corresponding iteration.
#'   \item \code{$orders} a multidimensional array where each slice is a matrix and represent an iteration. The row of each matrix correspond the order associated to the corresponding cluster.
#'   \item \code{$time} computational time.
#'   \item \code{$lkl} a matrix where each row is the likelihood of each observation computed at the corresponding iteration.
#'   \item \code{$norm_vec} a vector containing the normalization constant computed at the beginning of the algorithm.
#'   \item \code{I0_MCMC} traceplot for \eqn{I_0}.
#'   \item \code{I0_MCMC_01} a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{I_0} was accepted, \eqn{0} otherwise.
#'   \item \code{$kernel_ts} if TRUE data are time series.
#'   \item \code{$kernel_epi} if TRUE data are survival function.
#'   \item \code{$univariate_ts} TRUE if data is an univariate time series, FALSE if it is a multivariate time series.
#' }
#'
#'
#' @examples
#'
#'\donttest{
#' ## Univariate time series
#'
#' data_mat <- matrix(NA, nrow = 5, ncol = 100)
#'
#' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
#' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#' data_mat[4,] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
#' data_mat[5,] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
#'
#' out <- clust_cp(data = data_mat, n_iterations = 5000, n_burnin = 1000,
#'                  L = 1, params = list(phi = 0.5), B = 1000, kernel = "ts")
#'
#' print(out)
#'
#' ## Multivariate time series
#'
#'
#' data_array <- array(data = NA, dim = c(3,100,5))
#'
#' data_array[1,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_array[2,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_array[3,,1] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#'
#' data_array[1,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_array[2,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_array[3,,2] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#'
#' data_array[1,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#' data_array[2,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#' data_array[3,,3] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#'
#' data_array[1,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
#' data_array[2,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
#' data_array[3,,4] <- as.numeric(c(rnorm(25,0,0.135), rnorm(75,1,0.225)))
#'
#' data_array[1,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
#' data_array[2,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
#' data_array[3,,5] <- as.numeric(c(rnorm(25,0,0.155), rnorm(75,1,0.280)))
#'
#' out <- clust_cp(data = data_array, n_iterations = 3000, n_burnin = 1000,
#'                 params = list(phi = 0.5, k_0 = 0.25,
#'                               nu_0 = 5, S_0 = diag(0.1,3,3),
#'                               m_0 = rep(0,3)), B = 1000,  kernel = "ts")
#'
#' print(out)
#'
#' ## Epidemiological data
#'
#'
#'
#' data_mat <- matrix(NA, nrow = 5, ncol = 50)
#'
#' betas <- list(c(rep(0.45, 25),rep(0.14,25)),
#'               c(rep(0.55, 25),rep(0.11,25)),
#'               c(rep(0.50, 25),rep(0.12,25)),
#'               c(rep(0.52, 10),rep(0.15,40)),
#'               c(rep(0.53, 10),rep(0.13,40)))
#'
#' inf_times <- list()
#'
#' for(i in 1:5){
#'
#'   inf_times[[i]] <- sim_epi_data(10000, 10, 50, betas[[i]], 1/8)
#'
#'   vec <- rep(0,50)
#'   names(vec) <- as.character(1:50)
#'
#'   for(j in 1:50){
#'     if(as.character(j) %in% names(table(floor(inf_times[[i]])))){
#'       vec[j] = table(floor(inf_times[[i]]))[which(names(table(floor(inf_times[[i]]))) == j)]
#'     }
#'   }
#'   data_mat[i,] <- vec
#' }
#'
#' out <- clust_cp(data = data_mat, n_iterations = 100, n_burnin = 10,
#'                 params = list(M = 100, xi = 1/8), B = 1000, kernel = "epi")
#'
#' print(out)
#' }
#'
#' @references
#'
#' Corradin, R., Danese, L., KhudaBukhsh, W. R., & Ongaro, A. (2024). \emph{Model-based clustering of time-dependent observations with common structural changes}. arXiv preprint arXiv:2410.09552.
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
      if((!is.null(params$params) && !is.list(params))) stop("params must be a list")


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

      # save output

      result <- ClustCpObj(data = data_input,
                           n_iterations = n_iterations_input,
                           n_burnin = n_burnin_input,
                           clust = out$clust,
                           orders = out$orders,
                           time = out$time,
                           lkl = out$lkl,
                           norm_vec = out$norm_vec,
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
      phi_input = ifelse(is.null(params$phi), 1, params$phi)
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

      result <- ClustCpObj(data = data_input,
                           n_iterations = n_iterations_input,
                           n_burnin = n_burnin_input,
                           clust = out$clust,
                           orders = out$orders,
                           time = out$time,
                           lkl = out$lkl,
                           norm_vec = out$norm_vec,
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

    result <- ClustCpObj(data = data_input,
                         n_iterations = n_iterations_input,
                         n_burnin = n_burnin_input,
                         clust = out$clust,
                         orders = out$orders,
                         time = out$time,
                         lkl = out$lkl,
                         norm_vec = out$norm_vec,
                         I0_MCMC = out$I0_MCMC,
                         I0_MCMC_01 = out$I0_MCMC_01,
                         kernel_ts = FALSE,
                         kernel_epi = TRUE,
                         univariate_ts = TRUE)

  }

  return(result)


}
