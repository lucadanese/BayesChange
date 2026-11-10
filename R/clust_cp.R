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
#' @param user_seed seed for random number generation (default 1234). The output is fully determined by \code{user_seed}: two calls with the same data, arguments and \code{user_seed} return identical results, regardless of any previous call to \code{set.seed()}. The state of R's random number generator is restored on exit. If \code{NULL}, the current state of R's random number generator is used.
#' @param kernel can be "ts" if data are time series or "epi" if data are epidemic diffusions.
#'
#' @param params the model parameters, created with [params_uni()] (univariate
#'   time series), [params_multi()] (multivariate time series) or [params_epi()]
#'   (epidemic diffusions). Parameters not given take their default values, see
#'   [params_uni()]. A plain list with the same names is also accepted.
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
#' out <- clust_cp(data = stock_uni[1:5,], n_iterations = 2000, n_burnin = 500,
#'                 L = 1, q = 0.5, B = 1000,
#'                 params = params_uni(c = 1, phi = 0.1), kernel = "ts")
#'
#' print(out)
#'
#' ## Multivariate time series
#'
#' data("stock_multi")
#'
#' out <- clust_cp(data = stock_multi[,,1:5], n_iterations = 2000, n_burnin = 500,
#'                 L = 1, B = 1000,
#'                 params = params_multi(k_0 = 1, nu_0 = 10, S_0 = diag(1, 2, 2)),
#'                 kernel = "ts")
#'
#' print(out)
#'
#' ## Epidemic diffusions
#'
#' data("epi_synthetic_multi")
#'
#' out <- clust_cp(epi_synthetic_multi, n_iterations = 2000, n_burnin = 500,
#'                 L = 1, B = 1000, alpha_SM = 1,
#'                 params = params_epi(M = 100, I0_var = 0.1), kernel = "epi")
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
  if((!is.null(q)) && ((q >= 1) | (q <= 0))) stop("q must be in (0,1)")
  if((is.null(kernel) | !(kernel %in% c("ts", "epi")))) stop("kernel must be 'ts' or 'epi'")
  if((!is.null(print_progress) && (print_progress != TRUE & print_progress != FALSE))) stop("print_progress must be TRUE/FALSE")
  if(!is.null(user_seed) && !is.numeric(user_seed)) stop("user_seed must be an integer")

  user_seed <- local_user_seed(user_seed)

  if(kernel == "ts" && is.matrix(data)){

    rp <- resolve_params(params, "uni")
    # for backward compatibility, `q` given in params is still used
    report_sampler_settings(rp$sampler, "clust_cp", used = "q")
    if(!is.null(rp$sampler$q)) q <- rp$sampler$q
    pr <- rp$values

    out <- clust_cp_uni(data = data,
                        n_iterations = n_iterations,
                        B = B, L = L,
                        phi = pr$phi, a = pr$a, b = pr$b, c = pr$c,
                        q = q,
                        alpha_SM = alpha_SM,
                        print_progress = print_progress, user_seed = user_seed)

    result <- ClustCpObj(data = data,
                         n_iterations = n_iterations,
                         n_burnin = n_burnin,
                         clust = out$clust,
                         orders = out$orders,
                         time = out$time,
                         norm_vec = out$norm_vec,
                         entropy_MCMC = coda::mcmc(out$entropy),
                         lkl_MCMC = coda::mcmc(out$lkl),
                         kernel_ts = TRUE,
                         kernel_epi = FALSE,
                         univariate_ts = TRUE)

  } else if(kernel == "ts"){

    rp <- resolve_params(params, "multi", d = dim(data)[1])
    report_sampler_settings(rp$sampler, "clust_cp")
    pr <- rp$values

    out <- clust_cp_multi(data = data, n_iterations = n_iterations,
                          B = B, L = L, phi = pr$phi,
                          k_0 = pr$k_0, nu_0 = pr$nu_0,
                          S_0 = pr$S_0, m_0 = pr$m_0, q = q,
                          alpha_SM = alpha_SM,
                          print_progress = print_progress,
                          user_seed = user_seed)

    result <- ClustCpObj(data = data,
                         n_iterations = n_iterations,
                         n_burnin = n_burnin,
                         clust = out$clust,
                         orders = out$orders,
                         time = out$time,
                         norm_vec = out$norm_vec,
                         entropy_MCMC = coda::mcmc(out$entropy),
                         lkl_MCMC = coda::mcmc(out$lkl),
                         kernel_ts = TRUE,
                         kernel_epi = FALSE,
                         univariate_ts = FALSE)

  } else {

    rp <- resolve_params(params, "epi")
    report_sampler_settings(rp$sampler, "clust_cp")
    pr <- rp$values

    out <- clust_cp_epi(data = data, n_iterations = n_iterations,
                        M = pr$M, B = B, L = L,
                        xi = pr$xi, alpha_SM = alpha_SM, q = q,
                        a0 = pr$a0, b0 = pr$b0, I0_var = pr$I0_var,
                        avg_blk = pr$avg_blk, print_progress = print_progress,
                        user_seed = user_seed)

    result <- ClustCpObj(data = data,
                         n_iterations = n_iterations,
                         n_burnin = n_burnin,
                         clust = out$clust,
                         orders = out$orders,
                         time = out$time,
                         norm_vec = out$norm_vec,
                         entropy_MCMC = coda::mcmc(out$entropy),
                         lkl_MCMC = coda::mcmc(out$llik),
                         I0_MCMC = coda::mcmc(out$I0_MCMC),
                         kernel_ts = FALSE,
                         kernel_epi = TRUE,
                         univariate_ts = TRUE)

  }

  return(result)

}
