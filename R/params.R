#' Model parameters for detect_cp() and clust_cp()
#'
#' @description These functions create the list of model parameters passed to
#' the argument \code{params} of [detect_cp()] and [clust_cp()]:
#' \itemize{
#'   \item \code{params_uni()} for univariate time series;
#'   \item \code{params_multi()} for multivariate time series;
#'   \item \code{params_epi()} for epidemic diffusions.
#' }
#' Every argument has a default value, so only the parameters that differ from
#' the defaults need to be given. Printing the result shows all the values and
#' marks those left at their default.
#'
#' The settings of the sampling algorithms (\code{q}, \code{B}, \code{L},
#' \code{alpha_SM}) are not part of the model parameters: they are arguments of
#' [detect_cp()] and [clust_cp()].
#'
#' A plain list with the same names is still accepted by [detect_cp()] and
#' [clust_cp()]; unknown names raise a warning.
#'
#' @param a,b,c parameters of the Normal-Gamma prior of the mean and variance of
#'   each block: the variance has an Inverse-Gamma(\code{a}, \code{b}) prior and,
#'   given the variance \eqn{\sigma^2}, the mean has a
#'   \eqn{N(0, \sigma^2 / c)} prior.
#' @param m_0,k_0,nu_0,S_0 parameters of the Normal-Inverse-Wishart prior of the
#'   mean and covariance matrix of each block. If \code{NULL} (default),
#'   \code{m_0 = rep(0, d)}, \code{nu_0 = d + 1} and \code{S_0 = diag(0.1, d)},
#'   where \eqn{d} is the dimension of the time series.
#' @param phi autocorrelation parameter \eqn{\phi \in [0, 1)} of the
#'   Ornstein-Uhlenbeck process. Used by [clust_cp()], where it is fixed;
#'   [detect_cp()] estimates it.
#' @param prior_var_phi standard deviation of the Gaussian random-walk proposal
#'   on \eqn{\mathrm{logit}(\phi)} used in the Metropolis-Hastings update of
#'   \eqn{\phi}. Used by [detect_cp()].
#' @param prior_delta_c,prior_delta_d parameters of the shifted Gamma prior of
#'   the parameter \eqn{\delta} of the prior on the orders. Used by
#'   [detect_cp()].
#' @param M number of Monte Carlo replications used to compute the likelihood of
#'   an epidemic diffusion.
#' @param xi recovery rate, in (0, 1), fixed for each population at each time.
#' @param a0,b0 parameters of the Gamma prior of the infection rates, used in the
#'   integrated likelihood of the epidemic diffusions.
#' @param I0_var variance, in (0, 1), of the Metropolis-Hastings proposal for the
#'   proportion of infected individuals at time 0.
#' @param avg_blk prior average number of change points of each order. Used by
#'   [clust_cp()].
#' @param x an object created by \code{params_uni()}, \code{params_multi()} or
#'   \code{params_epi()}.
#' @param ... further arguments, ignored.
#'
#' @return An object of class \code{BayesChangeParams}: a named list with the
#'   values of the parameters.
#'
#' @examples
#'
#' # univariate time series: change only the prior of the block variance
#' params_uni(a = 2, b = 0.5)
#'
#' # multivariate time series: m_0, nu_0 and S_0 are set from the data
#' params_multi(k_0 = 1)
#'
#' # epidemic diffusions
#' params_epi(M = 100, I0_var = 0.1)
#'
#' @name params_uni
#' @export
params_uni <- function(a = 1, b = 1, c = 0.1,
                       phi = 0.1,
                       prior_var_phi = 0.1,
                       prior_delta_c = 1, prior_delta_d = 1){
  values <- list(a = a, b = b, c = c, phi = phi, prior_var_phi = prior_var_phi,
                 prior_delta_c = prior_delta_c, prior_delta_d = prior_delta_d)
  check_positive(values, c("a", "b", "c", "prior_var_phi", "prior_delta_c", "prior_delta_d"))
  check_phi(phi)
  new_params(values, "uni", given = names(match.call())[-1])
}

#' @rdname params_uni
#' @export
params_multi <- function(m_0 = NULL, k_0 = 0.5, nu_0 = NULL, S_0 = NULL,
                         phi = 0.1,
                         prior_var_phi = 0.1,
                         prior_delta_c = 1, prior_delta_d = 1){
  values <- list(m_0 = m_0, k_0 = k_0, nu_0 = nu_0, S_0 = S_0, phi = phi,
                 prior_var_phi = prior_var_phi,
                 prior_delta_c = prior_delta_c, prior_delta_d = prior_delta_d)
  check_positive(values, c("k_0", "nu_0", "prior_var_phi", "prior_delta_c", "prior_delta_d"))
  check_phi(phi)
  if(!is.null(m_0) && !(is.numeric(m_0) && is.null(dim(m_0)))) stop("m_0 must be a numeric vector", call. = FALSE)
  if(!is.null(S_0) && !(is.matrix(S_0) && nrow(S_0) == ncol(S_0))) stop("S_0 must be a square matrix", call. = FALSE)
  if(!is.null(m_0) && !is.null(S_0) && length(m_0) != nrow(S_0)) stop("m_0 and S_0 must have the same dimension", call. = FALSE)
  new_params(values, "multi", given = names(match.call())[-1])
}

#' @rdname params_uni
#' @export
params_epi <- function(M = 250, xi = 1/8, a0 = 4, b0 = 10, I0_var = 0.01,
                       avg_blk = 2){
  values <- list(M = M, xi = xi, a0 = a0, b0 = b0, I0_var = I0_var, avg_blk = avg_blk)
  check_positive(values, c("a0", "b0"))
  check_unit_interval(values, c("xi", "I0_var"))
  if(!(is.numeric(M) && length(M) == 1 && M >= 1)) stop("M must be at least equal to 1", call. = FALSE)
  if(!(is.numeric(avg_blk) && length(avg_blk) == 1 && avg_blk >= 0)) stop("avg_blk must be non-negative", call. = FALSE)
  new_params(values, "epi", given = names(match.call())[-1])
}

#' @rdname params_uni
#' @export
print.BayesChangeParams <- function(x, ...){
  type <- attr(x, "type")
  defaults <- attr(x, "defaults")
  title <- switch(type,
                  uni = "univariate time series",
                  multi = "multivariate time series",
                  epi = "epidemic diffusions")
  cat("BayesChange model parameters (", title, ")\n", sep = "")
  data_dependent <- c(m_0 = "rep(0, d)", nu_0 = "d + 1", S_0 = "diag(0.1, d)")
  for(nm in names(x)){
    v <- x[[nm]]
    shown <- if(is.null(v)) {
      paste0(data_dependent[[nm]], ", with d the dimension of the data")
    } else if(is.matrix(v)) {
      paste0(nrow(v), " x ", ncol(v), " matrix")
    } else {
      paste(format(v), collapse = " ")
    }
    cat(sprintf("  %-14s %s%s\n", nm, shown, if(nm %in% defaults) "  (default)" else ""))
  }
  invisible(x)
}

# --- internal helpers ----------------------------------------------------------

new_params <- function(values, type, given){
  structure(values, class = "BayesChangeParams", type = type,
            defaults = setdiff(names(values), given))
}

check_positive <- function(values, names){
  for(nm in names){
    v <- values[[nm]]
    if(!is.null(v) && !(is.numeric(v) && length(v) == 1 && v > 0))
      stop(nm, " must be a positive number", call. = FALSE)
  }
}

check_unit_interval <- function(values, names){
  for(nm in names){
    v <- values[[nm]]
    if(!(is.numeric(v) && length(v) == 1 && v > 0 && v < 1))
      stop(nm, " must be in (0, 1)", call. = FALSE)
  }
}

check_phi <- function(phi){
  if(!(is.numeric(phi) && length(phi) == 1 && phi >= 0 && phi < 1))
    stop("phi must be in [0, 1)", call. = FALSE)
}

# Settings of the samplers that used to be accepted in `params`.
sampler_settings <- c("q", "alpha_SM", "B", "L")

# Turn `params` (a BayesChangeParams object or a plain list) into the complete
# list of model parameters of the given type, filling the defaults that depend
# on the dimension `d` of the data. Returns the list of values and the sampler
# settings found in a plain list (for backward compatibility).
resolve_params <- function(params, type, d = NULL){

  constructor <- switch(type, uni = params_uni, multi = params_multi, epi = params_epi)
  constructor_name <- paste0("params_", type)
  sampler <- list()

  if(is.null(params)) params <- list()

  if(inherits(params, "BayesChangeParams")){
    if(attr(params, "type") != type){
      stop("params were created with params_", attr(params, "type"), "(), but these data require ",
           constructor_name, "()", call. = FALSE)
    }
  } else {
    if(!is.list(params)) stop("params must be created with ", constructor_name, "() or be a list", call. = FALSE)
    if(length(params) > 0 && (is.null(names(params)) || any(names(params) == ""))) stop("all elements of params must be named", call. = FALSE)

    known <- names(formals(constructor))
    in_sampler <- intersect(names(params), sampler_settings)
    unknown <- setdiff(names(params), c(known, sampler_settings))

    if(length(unknown) > 0){
      hints <- vapply(unknown, function(u){
        m <- agrep(u, known, max.distance = 0.2, value = TRUE)
        if(length(m) > 0) paste0("\"", u, "\" (did you mean \"", m[1], "\"?)") else paste0("\"", u, "\"")
      }, character(1))
      warning("unknown parameter(s) ignored: ", paste(hints, collapse = ", "),
              ". See ?", constructor_name, " for the parameters of this model.", call. = FALSE)
    }

    if(length(in_sampler) > 0) sampler <- params[in_sampler]

    params <- do.call(constructor, params[intersect(names(params), known)])
  }

  values <- unclass(params)
  attr(values, "type") <- NULL
  attr(values, "defaults") <- NULL

  if(type == "multi"){
    if(is.null(values$m_0)) values$m_0 <- rep(0, d)
    if(is.null(values$nu_0)) values$nu_0 <- d + 1
    if(is.null(values$S_0)) values$S_0 <- diag(0.1, d, d)
    if(length(values$m_0) != d) stop("m_0 must have length ", d, ", the dimension of the time series", call. = FALSE)
    if(!all(dim(values$S_0) == c(d, d))) stop("S_0 must be a ", d, " x ", d, " matrix, the dimension of the time series", call. = FALSE)
  }

  list(values = values, sampler = sampler)
}

# Warn about sampler settings found in a plain `params` list. `used` are those
# still honoured for backward compatibility, the others are ignored.
report_sampler_settings <- function(sampler, fun, used = character()){
  if(length(sampler) == 0) return(invisible(NULL))
  quote_names <- function(x) paste0("`", x, "`", collapse = ", ")
  u <- intersect(names(sampler), used)
  ig <- setdiff(names(sampler), used)
  if(length(u) > 0)
    warning(quote_names(u), " in params is deprecated: pass it as an argument of ", fun, "().", call. = FALSE)
  if(length(ig) > 0)
    warning(quote_names(ig), " in params ignored: settings of the sampler are arguments of ", fun,
            "(), not model parameters.", call. = FALSE)
  invisible(NULL)
}
