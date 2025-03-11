#' @name clust_cp
#' @export clust_cp
#'
#' @title Clustering time dependent obsevations with common change points.
#' @description The \code{clust_cp} function cluster observations with common change points. Data can be time series or survival functions.
#'
#' @param data a vector or a matrix. If a vector the algorithm for
#' univariate time series is used. If a matrix, where rows are the observations
#' and columns are the times, then the algorithm for multivariate time series is used.
#'
#' @param n_iterations number of MCMC iterations.
#' @param n_burnin number of iterations that must be excluded when computing the posterior estimate.
#' @param print_progress If TRUE (default) print the progress bar.
#' @param user_seed seed for random distribution generation.
#' @param kernel can be "ts" if data are time series or "epi" if data are survival functions.
#'
#' @param params a list of parameters:
#'
#' If the time series is univariate the following must be specified:
#'
#' \itemize{
#'   \item \code{q} probability of a split in the split-merge proposal and acceleration step.
#'   \item \code{B} number of orders for the normalization constant.
#'   \item \code{L} number of split-merge steps for the proposal step.
#'   \item \code{alpha_SM} \eqn{\alpha} for the split-merge proposal and acceleration step.
#'   \item \code{gamma},\code{a},\code{b},\code{c} parameters of the integrated likelihood.
#' }
#'
#' If the time series is multivariate the following must be specified:
#'
#' \itemize{
#'   \item \code{q} probability of a split in the split-merge proposal and acceleration step.
#'   \item \code{B} number of orders for the normalization constant.
#'   \item \code{L} number of split-merge steps for the proposal step.
#'   \item \code{gamma},\code{k_0},\code{nu_0},\code{phi_0},\code{m_0} parameters of the integrated likelihood.
#' }
#'
#' If data are survival functions:
#'
#' \itemize{
#'   \item \code{q} probability of a split in the split-merge proposal and acceleration step.
#'   \item \code{B} number of orders for the normalization constant.
#'   \item \code{L} number of split-merge steps for the proposal step.
#'   \item \code{alpha_SM} \eqn{\alpha} for the split-merge proposal and acceleration step.
#'   \item \code{M} number of Monte Carlo iterations when computing the likelihood of the survival function.
#'   \item \code{gamma} recovery rate fixed constant for each population at each time.
#'   \item \code{alpha} \eqn{\alpha} for the acceptance ration in the split-merge procedure.
#'   \item \code{dt},\code{a0},\code{b0},\code{c0},\code{d0} parameters for the computation of the integrated likelihood of the survival functions.
#'   \item \code{MH_var} variance for the Metropolis-Hastings estimation of the proportion of infected at time 0.
#'   \item \code{S0},\code{R0} parameters for the SDE solver.
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
#'   \item \code{$norm_vec} a vector containing the normalisation constant computed at the beginning of the algorithm.
#'   \item \code{$rho} a vector with the final estimate of the proportion of infected individuals at time 0.
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
#'                 params = list(L = 1, B = 1000, gamma = 0.5), kernel = "ts")
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
#'                 params = list(B = 1000, L = 1, gamma = 0.5, k_0 = 0.25,
#'                               nu_0 = 5, phi_0 = diag(0.1,3,3),
#'                               m_0 = rep(0,3)), kernel = "ts")
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
#'                 params = list(M = 100, L = 1, B = 1000), kernel = "epi")
#'
#' print(out)
#' }
#'
#' @references
#'
#' Corradin, R., Danese, L., KhudaBukhsh, W. R., & Ongaro, A. (2024). \emph{Model-based clustering of time-dependent observations with common structural changes}. arXiv preprint arXiv:2410.09552.
#'
#'
#'
clust_cp <- function(data,
                     n_iterations,
                     n_burnin = 0,
                     params = list(),
                     print_progress = TRUE,
                     user_seed = 1234,
                     kernel){

  if((!is.matrix(data)) && (!is.array(data))) stop("data must be a matrix or an array")
  if(is.null(n_iterations)) stop("missing number of iterations")
  if(n_iterations < 1) stop("number of iterations must be positive and > 1")
  if((is.null(kernel) | !(kernel %in% c("ts", "epi")))) stop("kernel must be 'ts' or 'epi'")
  if(n_burnin < 0) stop("number of burn in iterations must be positive and > 1")

  if(kernel == "ts"){

    if(is.matrix(data)){

      if((!is.null(params$B)) && (params$B < 1)) stop("params$B must be at least equal to 1")
      if((!is.null(params$L)) && (params$L < 1)) stop("params$L must be at least equal to 1")
      if((!is.null(params$gamma)) && ((params$gamma > 1) | (params$gamma < 0))) stop("params$gamma must be in (0,1)")
      if((!is.null(params$a)) && (params$a <= 0)) stop("params$a must be positive")
      if((!is.null(params$b)) && (params$b <= 0)) stop("params$b must be positive")
      if((!is.null(params$c)) && (params$c <= 0)) stop("params$c must be positive")
      if((!is.null(params$q)) && ((params$q >= 1) | (params$q <= 0))) stop("params$q must be in (0,1)")
      if((!is.null(params$alpha_SM)) && (params$alpha_SM <= 0)) stop("params$alpha_SM must be positive")
      if((!is.null(params$params) && !is.list(params))) stop("params must be a list")
      if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
      if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

      # substitute missing parameters with default

      n_burnin_input = n_burnin
      B_input = ifelse(is.null(params$B), 1000, params$B)
      L_input = ifelse(is.null(params$L), 1, params$L)
      gamma_input = ifelse(is.null(params$gamma), 1, params$gamma)
      a_input = ifelse(is.null(params$a), 1, params$a)
      b_input = ifelse(is.null(params$b), 1, params$b)
      c_input = ifelse(is.null(params$c), 1, params$c)
      q_input = ifelse(is.null(params$q), 0.5, params$q)
      alpha_SM_input = ifelse(is.null(params$alpha_SM), 0.1, params$alpha_SM)
      print_progress_input = print_progress
      user_seed_input = user_seed
      #

      data_input = data
      n_iterations_input = n_iterations

      out <- clust_cp_uni(data = data_input,
                          n_iterations = n_iterations_input,
                          B = B_input, L = L_input,
                          gamma = gamma_input, a = a_input, b = b_input,
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
                           rho = out$rho,
                           kernel_ts = TRUE,
                           kernel_epi = FALSE,
                           univariate_ts = TRUE)

    } else if(is.array(data)){

      if((!is.null(params$B)) && (params$B < 1)) stop("params$B must be at least equal to 1")
      if((!is.null(params$L)) && (params$L < 1)) stop("params$L must be at least equal to 1")
      if((!is.null(params$gamma)) && ((params$gamma > 1) | (params$gamma < 0))) stop("params$gamma must be in (0,1)")
      if((!is.null(params$k_0)) && (params$k_0 < 0)) stop("params$k_0 must be positive")
      if((!is.null(params$nu_0)) && (params$nu_0 < 0)) stop("params$nu_0 must be positive")
      if((!is.null(params$phi_0)) && nrow(params$phi_0) != ncol(params$phi_0)) stop("number of rows and columns must be the same in params$phi_0")
      if((!is.null(params$m_0)) && (length(params$m_0) != nrow(data))) stop("number of elements in params$m_0 must equal the number of observations")
      if((!is.null(params$q)) && ((params$q >= 1) | (params$q <= 0))) stop("params$q must be in (0,1)")
      if((!is.null(params$alpha_SM)) && (params$alpha_SM <= 0)) stop("params$alpha_SM must be positive")
      if((!is.null(params$params) && !is.list(params))) stop("params must be a list")
      if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
      if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

      # substitute missing parameters with default

      n_burnin_input = n_burnin
      B_input = ifelse(is.null(params$B), 1000, params$B)
      L_input = ifelse(is.null(params$L), 1, params$L)
      gamma_input = ifelse(is.null(params$gamma), 1, params$gamma)
      k_0_input = ifelse(is.null(params$k_0), 0.5, params$k_0)
      nu_0_input = ifelse(is.null(params$nu_0), nrow(data)+1, params$nu_0)
      q_input = ifelse(is.null(params$q), 0.5, params$q)
      alpha_SM_input = ifelse(is.null(params$alpha_SM), 0.1, params$alpha_SM)
      print_progress_input = print_progress
      user_seed_input = user_seed


      # with object matrix ifelse does not work
      if(is.null(params$phi_0)){phi_0_input = diag(0.1, nrow(data), nrow(data))} else{phi_0_input = params$phi_0}
      if(is.null(params$m_0)){m_0_input = rep(0, nrow(data))} else{m_0_input = params$m_0}
      #


      #

      data_input = data
      n_iterations_input = n_iterations

      out <- clust_cp_multi(data = data_input, n_iterations = n_iterations_input,
                            B = B_input, L = L_input, gamma = gamma_input,
                            k_0 = k_0_input, nu_0 = nu_0_input,
                            phi_0 = phi_0_input, m_0 = m_0_input, q = q_input,
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
                           rho = out$rho,
                           kernel_ts = TRUE,
                           kernel_epi = FALSE,
                           univariate_ts = FALSE)

    }

  }

  if(kernel == "epi"){

    if((!is.null(params$M)) && (params$M < 1)) stop("params$M must be at least equal to 1")
    if((!is.null(params$B)) && (params$B < 1)) stop("params$B must be at least equal to 1")
    if((!is.null(params$L)) && (params$L < 1)) stop("params$L must be at least equal to 1")
    if((!is.null(params$gamma)) && ((params$gamma <= 0) | (params$gamma >= 1))) stop("params$gamma must be in (0,1)")
    if((!is.null(params$alpha)) && (params$alpha <= 0)) stop("params$alpha must be positive")
    if((!is.null(params$q)) && ((params$q <= 0) | (params$q >= 1))) stop("params$q must be in (0,1)")
    if((!is.null(params$dt)) && ((params$dt <= 0) | (params$dt >= 1))) stop("params$dt must be in (0,1)")
    if((!is.null(params$a0)) && (params$a0 <= 0)) stop("params$a0 must be positive")
    if((!is.null(params$b0)) && (params$b0 <= 0)) stop("params$b0 must be positive")
    if((!is.null(params$c0)) && (params$c0 <= 0)) stop("params$c0 must be positive")
    if((!is.null(params$c0)) && (params$d0 <= 0)) stop("params$d0 must be positive")
    if((!is.null(params$MH_var)) && ((params$MH_var <= 0) | (params$MH_var >= 1))) stop("params$MH_var must be in (0,1)")
    if((!is.null(params$S0)) && (params$S0 <= 0)) stop("params$S0 must be positive")
    if((!is.null(params$R0)) && (params$R0 < 0)) stop("params$R0 must be at least 0")
    if((!is.null(params$dt)) && ((params$p <= 0) | (params$p >= 1))) stop("params$p must be in (0,1)")
    if((!is.null(params) && !is.list(params))) stop("params must be a list")
    if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
    if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

    # substitute missing parameters with default

    n_burnin_input = ifelse(is.null(n_burnin), 0, n_burnin)
    M_input = ifelse(is.null(params$M), 250, params$M)
    B_input = ifelse(is.null(params$M), 100, params$M)
    L_input = ifelse(is.null(params$L), 1, params$L)
    gamma_input = ifelse(is.null(params$gamma), 1/8, params$gamma)
    alpha_input = ifelse(is.null(params$alpha), 1, params$alpha)
    q_input = ifelse(is.null(params$q), 0.5, params$q)
    dt_input = ifelse(is.null(params$dt), 0.1, params$dt)
    a0_input = ifelse(is.null(params$a0), 4, params$a0)
    b0_input = ifelse(is.null(params$b0), 10, params$b0)
    c0_input = ifelse(is.null(params$c0), 1, params$c0)
    d0_input = ifelse(is.null(params$d0), 1, params$d0)
    MH_var_input = ifelse(is.null(params$MH_var), 0.01, params$MH_var)
    S0_input = ifelse(is.null(params$S0), 1, params$S0)
    R0_input = ifelse(is.null(params$R0), 0, params$R0)
    p_input = ifelse(is.null(params$p), 0.003, params$p)
    print_progress_input = print_progress
    user_seed_input = user_seed

    #

    data_input = data
    n_iterations_input = n_iterations

    out <- clust_cp_epi(data = data_input, n_iterations = n_iterations_input,
                        M = M_input, B = B_input, L = L_input,
                        gamma = gamma_input,alpha = alpha_input, q = q_input,
                        dt = dt_input, a0 = a0_input, b0 = b0_input,
                        c0 =  c0_input, d0 = d0_input, MH_var = MH_var_input,
                        S0 = S0_input, R0 = R0_input, p = p_input,
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
                         rho = out$rho,
                         kernel_ts = FALSE,
                         kernel_epi = TRUE,
                         univariate_ts = TRUE)

  }

  return(result)


}
