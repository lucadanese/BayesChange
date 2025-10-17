#' ClustCpObj class constructor
#'
#' @description A constructor for the \code{ClustCpObj} class. The class \code{ClustCpObj} contains...
#'
#' @param data a vector or a matrix containing the values of the time series;
#' @param n_iterations number of iterations of the MCMC algorithm;
#' @param n_burnin number of MCMC iterations to exclude in the posterior estimate;
#' @param clust a matrix with the clustering of each iteration.
#' @param orders a matrix where each row corresponds to the output order of the corresponding iteration;
#' @param time computational time in seconds;
#' @param lkl a vector with the likelihood of the final partition
#' @param norm_vec a vector with the estimated normalization constant.
#' @param I0_MCMC traceplot for \eqn{I_0}.
#' @param I0_MCMC_01 a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{I_0} was accepted, \eqn{0} otherwise.
#' @param kernel_ts if TRUE data are time series.
#' @param kernel_epi if TRUE data are epidemic diffusions.
#' @param univariate_ts TRUE/FALSE if time series is univariate or not;
#'
#'
#' @export
#'
ClustCpObj <- function(data = NULL,
                       n_iterations = NULL,
                       n_burnin = NULL,
                       clust = NULL,
                       orders = NULL,
                       time = NULL,
                       lkl = NULL,
                       norm_vec = NULL,
                       I0_MCMC = NULL,
                       I0_MCMC_01 = NULL,
                       kernel_ts = NULL,
                       kernel_epi = NULL,
                       univariate_ts = NULL){

  value <- list(data = data,
                n_iterations = n_iterations,
                n_burnin = n_burnin,
                clust = clust,
                orders = orders,
                time = time,
                lkl = lkl,
                norm_vec = norm_vec,
                I0_MCMC = I0_MCMC,
                I0_MCMC_01 = I0_MCMC_01,
                kernel_ts = kernel_ts,
                kernel_epi = kernel_epi,
                univariate_ts = univariate_ts)
  attr(value, "class") <- "ClustCpObj"
  value
}


#' ClustCpObj print method
#'
#' @description The \code{ClustCpObj} method prints which algorithm was run.
#' @param x an object of class \code{ClustCpObj}.
#' @param ... parameter of the generic method.
#'
#' @examples
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
#'                 params = list(L = 1, B = 1000, phi = 0.5), kernel = "ts")
#'
#' print(out)
#'
#' @rdname print.ClustCpObj
#' @export
#'
print.ClustCpObj <- function(x, ...) {
  cat("ClustCpObj object\n")
  if(x$kernel_ts){
    if(x$univariate){
      cat("Type: clustering univariate time series with common change points")
    } else {
      cat("Type: clustering multivariate time series with common change points")
    }
  } else if(x$kernel_epi){
    cat("Type: clustering epidemic diffusions with common change points")
  }
}


#' ClustCpObj summary method
#'
#' @description The \code{ClustCpObj} method returns a summary of the algorithm.
#' @param object an object of class \code{ClustCpObj}.
#' @param ... parameter of the generic method.
#'
#' @examples
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
#'                 params = list(L = 1, B = 1000, phi = 0.5), kernel = "ts")
#'
#' summary(out)
#'
#' @rdname summary.ClustCpObj
#' @export
#'
summary.ClustCpObj <- function(object, ...) {
  cat("ClustCpObj object\n")
  if(object$kernel_ts){
    if(object$univariate){
      cat("Clustering univariate time series:\n",
          "Number of burn-in iterations:", object$n_burnin, "\n",
          "Number of MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
          "Computational time:", round(object$time, digits = 2), "seconds\n")
    } else {
      cat("Clustering ", paste0(nrow(object$data[,,1]),"-dimensional time series:\n"),
          "Number of burn-in iterations:", object$n_burnin, "\n",
          "Number of MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
          "Computational time:", round(object$time, digits = 2), "seconds\n")
    }
  } else if (object$kernel_epi){

    cat("Clustering epidemic diffusions:\n",
        "Number of burn-in iterations:", object$n_burnin, "\n",
        "Number of MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
        "Computational time:", round(object$time, digits = 2), "seconds\n")

  }
}



#' Estimate the change points of the data
#'
#' @description  The \code{posterior_estimate} method estimates the change points of the data making use of the salso algorithm, for a \code{DetectCPObj} class object.
#'
#' @param object an object of class \code{ClustCpObj}.
#' @param loss The loss function used to estimate the final partition, it can be "VI", "binder", "omARI", "NVI", "ID", "NID".
#' @param maxNClusters maximum number of clusters in salso procedure.
#' @param nRuns number of runs in salso procedure.
#' @param maxZealousAttempts maximum number of zealous attempts in salso procedure.
#' @param ... parameter of the generic method.
#'
#' @details
#'
#' put details here
#'
#' @return
#'
#' The function returns a vector with the cluster assignment of each observation.
#'
#' @references
#'
#' #' D. B. Dahl, D. J. Johnson, and P. Müller (2022), Search Algorithms and Loss
#' Functions for Bayesian Clustering, \emph{Journal of Computational and
#' Graphical Statistics}, 31(4), 1189-1201, \doi{10.1080/10618600.2022.2069779}.
#'
#' @examples
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
#'                 params = list(L = 1, B = 1000, phi = 0.5), kernel = "ts")
#'
#' posterior_estimate(out)
#'
#' @rdname posterior_estimate.ClustCpObj
#' @export
#'
posterior_estimate.ClustCpObj <- function(object,
                               loss = "VI",
                               maxNClusters = 0,
                               nRuns = 16,
                               maxZealousAttempts = 10, ...) {

  mcmc_chain <- object$clust[(object$n_burnin + 1):object$n_iterations,]

  if(loss == "VI"){

    est_cp <- salso::salso(mcmc_chain, loss = "VI",
                           maxNClusters = maxNClusters,
                           nRuns = nRuns,
                           maxZealousAttempts = maxZealousAttempts)

    output <- as.numeric(est_cp)

    return(output)

  } else if(loss == "binder"){

    est_cp <- salso::salso(mcmc_chain, loss = "binder",
                           maxNClusters = maxNClusters,
                           nRuns = nRuns,
                           maxZealousAttempts = maxZealousAttempts)

    output <- as.numeric(est_cp)

    return(output)
  } else if (loss == "omARI"){
    est_cp <- salso::salso(mcmc_chain, loss = "omARI",
                           maxNClusters = maxNClusters,
                           nRuns = nRuns,
                           maxZealousAttempts = maxZealousAttempts)

    output <- as.numeric(est_cp)

    return(output)

  } else if (loss == "NVI"){
    est_cp <- salso::salso(mcmc_chain, loss = "NVI",
                           maxNClusters = maxNClusters,
                           nRuns = nRuns,
                           maxZealousAttempts = maxZealousAttempts)

    output <- as.numeric(est_cp)

    return(output)

  } else if (loss == "ID"){
    est_cp <- salso::salso(mcmc_chain, loss = "ID",
                           maxNClusters = maxNClusters,
                           nRuns = nRuns,
                           maxZealousAttempts = maxZealousAttempts)

    output <- as.numeric(est_cp)

    return(output)

  } else if (loss == "NID"){
    est_cp <- salso::salso(mcmc_chain, loss = "NID",
                           maxNClusters = maxNClusters,
                           nRuns = nRuns,
                           maxZealousAttempts = maxZealousAttempts)

    output <- as.numeric(est_cp)

    return(output)

  }
}


#' Plot estimated partition
#'
#' @description  The \code{plot} method plots the estimates partition through the salso algorithm, for a \code{ClustCpObj} class object.
#'
#' @param x an object of class \code{ClustCpObj}.
#' @param loss The loss function used to estimate the final partition, it can be "VI", "binder", "omARI", "NVI", "ID", "NID".
#' @param maxNClusters maximum number of clusters in salso procedure.
#' @param nRuns number of runs in salso procedure.
#' @param maxZealousAttempts maximum number of zealous attempts in salso procedure.
#' @param y parameter of the generic method.
#' @param ... parameter of the generic method.
#'
#' @return
#'
#' The function returns a ggplot object representing the time series or the epidemic diffusions colored according to the final partition.
#'
#' @examples
#'
#'\donttest{
#' ## Time series
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
#'                 params = list(L = 1, B = 1000, phi = 0.5), kernel = "ts")
#'
#' plot(out)
#'
#'
#' ## Epidemic diffusions
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
#'   inf_times[[i]] <- sim_epi_data(10000, 10, 50, betas[[i]], 1/8)
#'   vec <- rep(0,50)
#'   names(vec) <- as.character(1:50)
#'   for(j in 1:50){
#'     if(as.character(j) %in% names(table(floor(inf_times[[i]])))){
#'       vec[j] = table(floor(inf_times[[i]]))[which(names(table(floor(inf_times[[i]]))) == j)]
#'     }
#'   }
#'   data_mat[i,] <- vec
#' }
#'
#' out <- clust_cp(data = data_mat, n_iterations = 100, n_burnin = 10,
#'                 params = list(M = 100, L = 1, B = 100), kernel = "epi")
#'
#' plot(out)
#' }
#'
#' @rdname plot.ClustCpObj
#' @export
#'
plot.ClustCpObj <- function(x, y = NULL,
                             loss = "VI",
                             maxNClusters = 0,
                             nRuns = 16,
                             maxZealousAttempts = 10, ...) {

  Observation <- Time <- Value <- Cluster <- Count <- y <- obs <- NULL

  if(x$kernel_ts){
    if(x$univariate_ts){

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      .data_plot <- x$data

      .data_plot <- as.data.frame(t(x$data))
      .data_plot <- tidyr::pivot_longer(.data_plot, cols = dplyr::everything(), names_to = "Observation", values_to = "Value")
      .data_plot <- dplyr::mutate(.data_plot, Observation = as.numeric(gsub("V", "", Observation)))
      .data_plot <- dplyr::group_by(.data_plot, Observation)
      .data_plot <- dplyr::mutate(.data_plot, Time = dplyr::row_number())
      .data_plot <-  dplyr::ungroup(.data_plot)

      .data_plot$Cluster <- rep(est_cp, max(.data_plot$Time))
      .data_plot$Cluster <- as.factor(.data_plot$Cluster)
      .data_plot$Observation <- as.factor(.data_plot$Observation)

      ggplot2::ggplot(.data_plot) +
        #ggplot2::geom_line(ggplot2::aes(x = Time, y = Value, color = Cluster, group = Observation, linetype = Observation)) +
        ggplot2::geom_line(ggplot2::aes(x = Time, y = Value, color = Observation, group = Observation, linetype = Cluster)) +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Value") +
        ggplot2::scale_colour_brewer(palette = "Set1") +
        ggplot2::theme_minimal()

    } else {

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)


      .data_plot <- data.frame(Value = numeric(0))

      count = 0
      for (i in 1:dim(x$data)[3]) {
        mat <- x$data[,,i]
        for(j in 1:nrow(mat)){
          count = count + 1
          .data_plot <- rbind(.data_plot, data.frame(Value = as.vector(mat[j,]),
                                                     Observation = i,
                                                     Cluster = est_cp[i],
                                                     Time = 1:ncol(mat),
                                                     Count = count))
        }
      }


      .data_plot$Observation <- as.factor(.data_plot$Observation)
      .data_plot$Cluster <- as.factor(.data_plot$Cluster)
      .data_plot$Count <- as.factor(.data_plot$Count)

      ggplot2::ggplot(.data_plot) +
        ggplot2::geom_line(ggplot2::aes(x = Time, y = Value, color = Observation, group = Count, linetype = Cluster)) +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Value") +
        ggplot2::scale_colour_brewer(palette = "Set1") +
        ggplot2::theme_minimal()

    }
  } else if(x$kernel_epi){

    est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

    .df_sf_plot <- data.frame(y = as.vector(sapply(1:nrow(x$data), function(y) 1 - cumsum(x$data[y,]) / sum(x$data[y,]))),
                             x = rep(1:ncol(x$data), nrow(x$data)),
                             obs = as.factor(rep(1:nrow(x$data), each = ncol(x$data))),
                             Cluster = as.factor(rep(est_cp, each = ncol(x$data))))


    ggplot2::ggplot(.df_sf_plot, ggplot2::aes(x = x, y = y, group = obs, colour = Cluster)) +
      ggplot2::geom_line(lwd = 0.5) +
      ggplot2::xlab("Time") +
      ggplot2::ylab("Proportion of Infected Individuals") +
      ggplot2::scale_colour_brewer(palette = "Set1") +
      ggplot2::theme_minimal()

  }

}
