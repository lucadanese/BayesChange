#' DetectCpObj class constructor
#'
#' @description A constructor for the \code{DetectCpObj} class. The class \code{DetectCpObj} contains...
#'
#' @param data a vector or a matrix containing the values of the time series;
#' @param n_iterations number of iterations of the MCMC algorithm;
#' @param n_burnin number of MCMC iterations to exclude in the posterior estimate;
#' @param orders a matrix where each row corresponds to the output order of the corresponding iteration;
#' @param time computational time in seconds;
#'
#'
#' @param gamma_MCMC traceplot for \eqn{\gamma};
#' @param gamma_MCMC_01 a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\gamma} was accepted, \eqn{0} otherwise;
#' @param sigma_MCMC traceplot for \eqn{\sigma};
#' @param sigma_MCMC_01 a \eqn{0/1} vector, the \eqn{n}-th element is equal to \eqn{1} if the proposed \eqn{\sigma} was accepted, \eqn{0} otherwise;
#' @param theta_MCMC traceplot for \eqn{\theta};
#' @param univariate_ts TRUE/FALSE if time series is univariate or not;
#'
#'
#' @export
#'
DetectCpObj <- function(data = NULL,
                         n_iterations = NULL,
                         n_burnin = NULL,
                         orders = NULL,
                         time = NULL,
                         gamma_MCMC = NULL,
                         gamma_MCMC_01 = NULL,
                         sigma_MCMC = NULL,
                         sigma_MCMC_01 = NULL,
                         theta_MCMC = NULL,
                         univariate_ts = NULL){

  value <- list(data = data,
                n_iterations = n_iterations,
                n_burnin = n_burnin,
                orders = orders,
                time = time,
                gamma_MCMC = gamma_MCMC,
                gamma_MCMC_01 = gamma_MCMC_01,
                sigma_MCMC = sigma_MCMC,
                sigma_MCMC_01 = sigma_MCMC_01,
                theta_MCMC = theta_MCMC,
                univariate_ts = univariate_ts)
  attr(value, "class") <- "DetectCpObj"
  value
}

#' DetectCpObj print method
#'
#' @description The \code{DetectCpObj} method prints which algorithm was run.
#' @param x an object of class \code{DetectCpObj}.
#' @param ... parameter of the generic method.
#'
#' @examples
#' data_mat <- matrix(NA, nrow = 3, ncol = 100)
#'
#' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
#' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#'
#' out <- detect_cp(data = data_mat, n_iterations = 2500, n_burnin = 500,
#'                 params = list(q = 0.25, k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
#'                               par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1))
#' print(out)
#'
#' @rdname print.DetectCpObj
#' @export
#'
print.DetectCpObj <- function(x, ...) {
  cat("DetectCpObj object\n")
  if(x$univariate){
    cat("Type: change points detection on univariate time series")
  } else {
    cat("Type: change points detection on multivariate time series")
  }
}

#' DetectCpObj summary method
#'
#' @description The \code{DetectCpObj} method returns a summary of the algorithm.
#' @param object an object of class \code{DetectCpObj};
#' @param ... parameter of the generic method.
#'
#' @examples
#'
#' data_mat <- matrix(NA, nrow = 3, ncol = 100)
#'
#' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
#' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#'
#' out <- detect_cp(data = data_mat, n_iterations = 2500, n_burnin = 500,
#'                 params = list(q = 0.25, k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3), m_0 = rep(0,3),
#'                               par_theta_c = 2, par_theta_d = 0.2, prior_var_gamma = 0.1))
#' summary(out)
#'
#' @rdname summary.DetectCpObj
#' @export
#'
summary.DetectCpObj <- function(object, ...) {
  cat("DetectCpObj object\n")
  if(object$univariate){
    cat("Detecting change points on an univariate time series:\n",
        "Number of burn-in iterations:", object$n_burnin, "\n",
        "Number of MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
        "Computational time:", round(object$time, digits = 2), "seconds\n")
  } else {
    cat("Detecting change points on a", paste0(nrow(object$data),"-dimensional time series:\n"),
        "Number of burn-in iterations:", object$n_burnin, "\n",
        "Number of MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
        "Computational time:", round(object$time, digits = 2), "seconds\n")
  }
}


#' set generic
#' @name posterior_estimate
#' @keywords internal
#' @export
#'
posterior_estimate <- function (object, ...) {
  UseMethod("posterior_estimate")
}

#' Estimate the change points of the data
#'
#' @description  The \code{posterior_estimate} method estimates the change points of the data making use of the salso algorithm, for a \code{DetectCPObj} class object.
#'
#' @param object an object of class \code{DetectCPObj}.
#' @param loss The loss function used to estimate the final partition, it can be "VI", "binder", "omARI", "NVI", "ID", "NID".
#' @param maxNClusters maximum number of clusters in salso procedure.
#' @param nRuns number of runs in salso procedure.
#' @param maxZealousAttempts maximum number of zealous attempts in salso procedure.
#' @param ... parameter of the generic method.
#'
#'
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
#' D. B. Dahl, D. J. Johnson, and P. Müller (2022), Search Algorithms and Loss
#' Functions for Bayesian Clustering, \emph{Journal of Computational and
#' Graphical Statistics}, 31(4), 1189-1201, \doi{10.1080/10618600.2022.2069779}.
#'
#' @examples
#'
#' data_vec <- as.numeric(c(rnorm(50,0,0.1), rnorm(50,1,0.25)))
#'
#'
#' out <- detect_cp(data = data_vec, n_iterations = 2500, n_burnin = 500,
#'                  params = list(q = 0.25, phi = 0.1, a = 1, b = 1, c = 0.1))
#'
#' posterior_estimate(out)
#'
#' @rdname posterior_estimate.DetectCpObj
#' @export
#'
posterior_estimate.DetectCpObj <- function(object,
                               loss = "VI",
                               maxNClusters = 0,
                               nRuns = 16,
                               maxZealousAttempts = 10,...) {

  mcmc_chain <- object$orders[(object$n_burnin + 1):object$n_iterations,]

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



#' Plot estimated change points
#'
#' @description  The \code{plot} method plots the estimates change points estimated through the salso algorithm, for a \code{DetectCpObj} class object.
#'
#' @param x an object of class \code{DetectCPObj}.
#' @param plot_freq if TRUE also the histogram with the empirical frequency of each change point is plotted.
#' @param loss The loss function used to estimate the final partition, it can be "VI", "binder", "omARI", "NVI", "ID", "NID".
#' @param maxNClusters maximum number of clusters in salso procedure.
#' @param nRuns number of runs in salso procedure.
#' @param maxZealousAttempts maximum number of zealous attempts in salso procedure.
#' @param y,... parameters of the generic method.
#'
#'
#'
#' @return
#'
#' The function returns a ggplot object representing the detected change points. If \code{plot_freq = TRUE} is plotted also an histogram with the frequency of times that a change point has been detected in the MCMC chain.
#'
#'
#' @examples
#'
#' data_mat <- matrix(NA, nrow = 3, ncol = 100)
#'
#' data_mat[1,] <- as.numeric(c(rnorm(50,0,0.100), rnorm(50,1,0.250)))
#' data_mat[2,] <- as.numeric(c(rnorm(50,0,0.125), rnorm(50,1,0.225)))
#' data_mat[3,] <- as.numeric(c(rnorm(50,0,0.175), rnorm(50,1,0.280)))
#'
#' out <- detect_cp(data = data_mat, n_iterations = 2500, n_burnin = 500,
#'                  params = list(q = 0.25, k_0 = 0.25, nu_0 = 4, phi_0 = diag(1,3,3),
#'                                m_0 = rep(0,3), par_theta_c = 2, par_theta_d = 0.2,
#'                                prior_var_gamma = 0.1))
#' plot(out)
#'
#'
#'
#' @rdname plot.DetectCpObj
#' @export
#'
plot.DetectCpObj <- function(x, y = NULL,
                             plot_freq = FALSE,
                             loss = "VI",
                             maxNClusters = 0,
                             nRuns = 16,
                             maxZealousAttempts = 10, ...) {


  time <- V2 <- y <- NULL

  if(plot_freq == FALSE){


    if(x$univariate_ts){

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      cp <- cumsum(table(est_cp))[-length(table(est_cp))]

      vec_data <- x$data

      .data_plot <- as.data.frame(cbind(vec_data))
      .data_plot$time <- rep(1:length(x$data))

      p1 <- ggplot2::ggplot(.data_plot) +
        ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data),  linetype = 1) +
        ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
        ggplot2::labs(x = "time",
                      y = " ",
                      color = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="top", legend.key.width = ggplot2::unit(1, 'cm'))

      p1


    } else {

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      cp <- cumsum(table(est_cp))[-length(table(est_cp))]

      vec_data <- as.numeric()

      for(i in 1:nrow(x$data)){
        vec_data <- c(vec_data,x$data[i,])
      }

      .data_plot <- as.data.frame(cbind(vec_data, sort(rep(1:nrow(x$data),ncol(x$data)))))
      .data_plot$V2 <- factor(.data_plot$V2, labels = unique(paste0("obs ", .data_plot$V2)) )
      .data_plot$time <- rep(1:ncol(x$data),nrow(x$data))

      p1 <- ggplot2::ggplot(.data_plot) +
        ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data, color = V2),  linetype = 1) +
        ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
        ggplot2::labs(x = "time",
                      y = " ",
                      color = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="top", legend.key.width = ggplot2::unit(1, 'cm'))

      p1

    }

  } else {


    if(x$univariate_ts){

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      cp <- cumsum(table(est_cp))[-length(table(est_cp))]

      vec_data <- x$data

      .data_plot <- as.data.frame(cbind(vec_data))
      .data_plot$time <- rep(1:length(x$data))

      p1 <- ggplot2::ggplot(.data_plot) +
        ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data),  linetype = 1) +
        ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
        ggplot2::labs(x = "time",
                      y = " ",
                      color = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="top", legend.key.width = ggplot2::unit(1, 'cm'))

      p1

      x_unique <- unique(.data_plot$time)
      b <- rep(0, length(x$data))

      for(i in 1:x$n_iterations){

        cp_iteration <- cumsum(table(x$orders[i,]))[-length(table(x$orders[i,]))]

        b[cp_iteration] = b[cp_iteration] + 1

      }

      b <- b/10000

      p2 <- ggplot2::ggplot(data.frame(x = x_unique, y =b)) +
        ggplot2::geom_bar(ggplot2::aes(x = x_unique, y = y), stat="identity", width = 0.5, col = "black") +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(angle = 90)) +
        ggplot2::scale_y_continuous(breaks = c(0,.5,1)) +
        ggplot2::ylab("Prob.") +
        ggplot2::xlab("Date") +
        ggplot2::theme_minimal()

      ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(2,1))

    } else {

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      cp <- cumsum(table(est_cp))[-length(table(est_cp))]

      vec_data <- as.numeric()

      for(i in 1:nrow(x$data)){
        vec_data <- c(vec_data,x$data[i,])
      }

      .data_plot <- as.data.frame(cbind(vec_data, sort(rep(1:nrow(x$data),ncol(x$data)))))
      .data_plot$V2 <- factor(.data_plot$V2, labels = unique(paste0("obs ", .data_plot$V2)) )
      .data_plot$time <- rep(1:ncol(x$data),nrow(x$data))

      p1 <- ggplot2::ggplot(.data_plot) +
        ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data, color = V2),  linetype = 1) +
        ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
        ggplot2::labs(x = "time",
                      y = " ",
                      color = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="top", legend.key.width = ggplot2::unit(1, 'cm'))


      x_unique <- unique(.data_plot$time)
      b <- rep(0, ncol(x$data))

      for(i in 1:x$n_iterations){

        cp_iteration <- cumsum(table(x$orders[i,]))[-length(table(x$orders[i,]))]

        b[cp_iteration] = b[cp_iteration] + 1

      }

      b <- b/10000

      p2 <- ggplot2::ggplot(data.frame(x = x_unique, y =b)) +
        ggplot2::geom_bar(ggplot2::aes(x = x_unique, y = y), stat="identity", width = 0.5, col = "black") +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(angle = 90)) +
        ggplot2::scale_y_continuous(breaks = c(0,.5,1)) +
        ggplot2::ylab("Prob.") +
        ggplot2::xlab("Date") +
        ggplot2::theme_minimal()

      ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(2,1))

    }

  }

}
