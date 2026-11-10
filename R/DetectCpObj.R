#' DetectCpObj class constructor
#'
#' @description
#' Constructor for the \code{DetectCpObj} class. This class stores the output
#' of the Bayesian change–point detection algorithm, including MCMC traces,
#' allocation orders, and computational information.
#'
#' @param data A vector or matrix containing the observed time series.
#' @param n_iterations Total number of MCMC iterations.
#' @param n_burnin Number of burn-in iterations to discard.
#' @param orders A matrix where each row corresponds to the latent block
#'        assignment (order) of the time indices at each MCMC iteration.
#' @param time Computational time in seconds.
#'
#' @param entropy_MCMC A \code{coda::mcmc} object containing MCMC samples of the entropy measure.
#' @param lkl_MCMC A \code{coda::mcmc} object containing MCMC samples of the log-likelihood.
#' @param phi_MCMC A \code{coda::mcmc} object containing MCMC draws for \eqn{\gamma}.
#' @param sigma_MCMC A \code{coda::mcmc} object containing MCMC draws for \eqn{\sigma}.
#' @param delta_MCMC A \code{coda::mcmc} object containing MCMC draws for \eqn{\delta}.
#' @param I0_MCMC A \code{coda::mcmc} object containing MCMC draws for \eqn{I_0}.
#'
#' @param kernel_ts Logical; TRUE if the model for time series data is used.
#' @param kernel_epi Logical; TRUE if the epidemic diffusion model is used.
#' @param univariate_ts Logical; TRUE if the time series is univariate, FALSE otherwise.
#'
#' @export
#'
DetectCpObj <- function(data = NULL,
                         n_iterations = NULL,
                         n_burnin = NULL,
                         orders = NULL,
                         time = NULL,
                         entropy_MCMC = NULL,
                         lkl_MCMC = NULL,
                         phi_MCMC = NULL,
                         sigma_MCMC = NULL,
                         delta_MCMC = NULL,
                         I0_MCMC = NULL,
                         kernel_ts = NULL,
                         kernel_epi = NULL,
                         univariate_ts = NULL){

  value <- list(data = data,
                n_iterations = n_iterations,
                n_burnin = n_burnin,
                orders = orders,
                time = time,
                entropy_MCMC = entropy_MCMC,
                lkl_MCMC = lkl_MCMC,
                phi_MCMC = phi_MCMC,
                sigma_MCMC = sigma_MCMC,
                delta_MCMC = delta_MCMC,
                I0_MCMC = I0_MCMC,
                kernel_ts = kernel_ts,
                kernel_epi = kernel_epi,
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
#'
#' data("eu_inflation")
#'
#' model_params <- params_uni(a = 1, b = 1, c = 1, prior_var_phi = 0.1,
#'                    prior_delta_c = 1, prior_delta_d = 1)
#'
#' out <- detect_cp(data = eu_inflation[1,], n_iterations = 1000,
#'                  n_burnin = 100, q = 0.5, params = model_params,
#'                  kernel = "ts")
#' print(out)
#'
#' @rdname print.DetectCpObj
#' @export
#'
print.DetectCpObj <- function(x, ...) {
  cat("DetectCpObj object\n")
  if(x$kernel_ts){
    if(x$univariate_ts){
      cat("Type: change points detection on univariate time series\n")
    } else {
      cat("Type: change points detection on multivariate time series\n")
    }
  }
  if(x$kernel_epi){
    cat("Type: change points detection on an epidemic diffusion\n")
  }

  invisible(x)

}

#' DetectCpObj summary method
#'
#' @description The \code{DetectCpObj} method returns a summary of the algorithm.
#' @param object An object of class \code{DetectCpObj};
#' @param ... parameter of the generic method.
#'
#' @examples
#'
#' data("eu_inflation")
#'
#' model_params <- params_uni(a = 1, b = 1, c = 1, prior_var_phi = 0.1,
#'                    prior_delta_c = 1, prior_delta_d = 1)
#'
#' out <- detect_cp(data = eu_inflation[1,], n_iterations = 1000,
#'                  n_burnin = 100, q = 0.5, params = model_params,
#'                  kernel = "ts")
#' summary(out)
#'
#' @rdname summary.DetectCpObj
#' @export
#'
summary.DetectCpObj <- function(object, ...) {
  cat("DetectCpObj object\n")
  if(object$kernel_ts){
    if(object$univariate){
      cat("Change point detection summary:\n",
          "- Data: univariate time series\n",
          "- Burn-in iterations:", object$n_burnin, "\n",
          "- MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
          "- Average number of detected change points:", round(mean(apply(object$orders[(object$n_burnin+1):object$n_iterations,], 1,  function(x) max(x))),2), "\n",
          "- Computational time:", round(object$time, 2), "seconds\n",
          "\nUse plot() for a detailed visualization or posterior_estimate() to analyze the detected change points.\n")
    } else {
      cat("Change point detection summary:\n",
          "- Data:", nrow(object$data), "-dimensional time series\n",
          "- Burn-in iterations:", object$n_burnin, "\n",
          "- MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
          "- Average number of detected change points:",
          round(mean(apply(object$orders[(object$n_burnin + 1):object$n_iterations, ], 1, function(x) max(x))), 2), "\n",
          "- Computational time:", round(object$time, 2), "seconds\n",
          "\nUse plot() for a detailed visualization or posterior_estimate() to analyze the detected change points.\n")

    }
  }

  if(object$kernel_epi){
    cat("Change point detection summary:\n",
        "- Data: epidemic diffusion\n",
        "- Burn-in iterations:", object$n_burnin, "\n",
        "- MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
        "- Average number of detected change points:",
        round(mean(apply(object$orders[(object$n_burnin + 1):object$n_iterations, ], 1, function(x) max(x))), 2), "\n",
        "- Computational time:", round(object$time, 2), "seconds\n",
        "\nUse plot() for a detailed visualization or posterior_estimate() to analyze the detected change points.\n")

  }


}


#' Point estimate of change points or of the partition
#'
#' @description \code{posterior_estimate} summarises the MCMC output of
#' [detect_cp()] or [clust_cp()] with a single point estimate, obtained by
#' minimising the posterior expected loss with the SALSO algorithm of
#' Dahl, Johnson and Müller (2022), as implemented in [salso::salso()].
#'
#' \itemize{
#'   \item For a \code{DetectCpObj} (output of [detect_cp()]) the sampled
#'   orders are treated as partitions of the time indices and the estimate is
#'   a partition of the times into contiguous blocks. With
#'   \code{show_cp = TRUE} the positions of the change points are returned
#'   instead.
#'   \item For a \code{ClustCpObj} (output of [clust_cp()]) the estimate is
#'   the partition of the observations (time series or epidemic diffusions)
#'   into clusters.
#' }
#'
#' @param object an object of class \code{DetectCpObj} or \code{ClustCpObj}.
#' @param show_cp (\code{DetectCpObj} only) logical; if \code{TRUE} return the
#'   positions of the estimated change points, i.e. the first time index of
#'   every block after the first one. If \code{FALSE} (default) return the
#'   block membership of each time index.
#' @param loss the loss function minimised by SALSO. One of
#'   \describe{
#'     \item{\code{"VI"}}{variation of information (default; Meilă, 2007);}
#'     \item{\code{"binder"}}{Binder's loss (Binder, 1978);}
#'     \item{\code{"omARI"}}{one minus the adjusted Rand index;}
#'     \item{\code{"NVI"}}{normalised variation of information;}
#'     \item{\code{"ID"}}{information distance;}
#'     \item{\code{"NID"}}{normalised information distance.}
#'   }
#'   See [salso::salso()] for details.
#' @param maxNClusters maximum number of blocks (\code{DetectCpObj}) or
#'   clusters (\code{ClustCpObj}) considered by SALSO. The default \code{0}
#'   sets no upper bound.
#' @param nRuns number of independent SALSO runs.
#' @param maxZealousAttempts maximum number of zealous updates per SALSO run.
#'   The default (\code{100}) is larger than the one of \pkg{salso}
#'   (\code{10}), because a time series can easily contain more than ten
#'   blocks; with too few attempts SALSO warns that \emph{"the number of
#'   possible zealous attempts exceeded the maximum"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return For a \code{DetectCpObj}, an integer vector with the block of each
#'   time index (blocks are numbered from 1 in temporal order) or, if
#'   \code{show_cp = TRUE}, the positions of the change points. For a
#'   \code{ClustCpObj}, an integer vector with the cluster of each observation.
#'
#' @references
#'
#' D. B. Dahl, D. J. Johnson, and P. Müller (2022), Search Algorithms and Loss
#' Functions for Bayesian Clustering, \emph{Journal of Computational and
#' Graphical Statistics}, 31(4), 1189-1201, \doi{10.1080/10618600.2022.2069779}.
#'
#' @examples
#'
#' data("eu_inflation")
#'
#' model_params <- params_uni(a = 1, b = 1, c = 1, prior_var_phi = 0.1,
#'                    prior_delta_c = 1, prior_delta_d = 1)
#'
#' out <- detect_cp(data = eu_inflation[1,], n_iterations = 1000,
#'                  n_burnin = 100, q = 0.5, params = model_params,
#'                  kernel = "ts")
#'
#' posterior_estimate(out, loss = "binder")
#' posterior_estimate(out, loss = "binder", show_cp = TRUE)
#'
#' data("stock_uni")
#'
#' model_params <- params_uni(a = 1, b = 1, c = 1, phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:3,], n_iterations = 1000, n_burnin = 100,
#'                 L = 1, q = 0.5, B = 500, params = model_params, kernel = "ts")
#'
#' posterior_estimate(out, loss = "binder")
#'
#' @export
posterior_estimate <- function (object, ...) {
  UseMethod("posterior_estimate")
}

# Run SALSO on a matrix of sampled partitions (one per row).
run_salso <- function(mcmc_chain, loss, maxNClusters, nRuns, maxZealousAttempts){
  loss <- match.arg(loss, c("VI", "binder", "omARI", "NVI", "ID", "NID"))
  # maxNClusters = 0: no upper bound (salso's own default caps the estimate at
  # the largest number of clusters in the sample and warns when it is reached)
  if(maxNClusters == 0) maxNClusters <- ncol(mcmc_chain)
  as.integer(salso::salso(mcmc_chain, loss = loss,
                          maxNClusters = maxNClusters,
                          nRuns = nRuns,
                          maxZealousAttempts = maxZealousAttempts))
}

#' @rdname posterior_estimate
#' @export
posterior_estimate.DetectCpObj <- function(object,
                               show_cp = FALSE,
                               loss = "VI",
                               maxNClusters = 0,
                               nRuns = 16,
                               maxZealousAttempts = 100, ...) {

  mcmc_chain <- object$orders[(object$n_burnin + 1):object$n_iterations, , drop = FALSE]

  est <- run_salso(mcmc_chain, loss, maxNClusters, nRuns, maxZealousAttempts)

  # SALSO is not aware of the temporal ordering, so the labels it returns need
  # not form contiguous blocks. Every maximal run of equal labels is a block.
  new_block <- c(TRUE, est[-1] != est[-length(est)])
  output <- cumsum(new_block)

  if(show_cp){
    output <- which(new_block)[-1]
  }

  return(output)
}

#' Plot estimated change points
#'
#' @description  The \code{plot} method plots the estimated change points estimated through the salso algorithm, for a \code{DetectCpObj} class object.
#'
#' @param x An object of class \code{DetectCpObj}.
#' @param plot_freq Logical; TRUE also the histogram with the empirical frequency of each change point is plotted.
#' @param loss The loss function used to estimate the final partition, it can be "VI", "binder", "omARI", "NVI", "ID", "NID".
#' @param maxNClusters Maximum number of clusters in salso procedure.
#' @param nRuns Number of runs in salso procedure.
#' @param maxZealousAttempts Maximum number of zealous attempts in salso procedure.
#' @param y,... parameters of the generic method.
#'
#'
#' @return
#'
#' The function returns a ggplot object representing the detected change points. If \code{plot_freq = TRUE} is plotted also an histogram with the frequency of times that a change point has been detected in the MCMC chain.
#'
#' @examples
#'
#' ## Univariate time series
#'
#' data("eu_inflation")
#'
#' model_params <- params_uni(a = 1, b = 1, c = 1, prior_var_phi = 0.1,
#'                    prior_delta_c = 1, prior_delta_d = 1)
#'
#' out <- detect_cp(data = eu_inflation[1,], n_iterations = 1000,
#'                  n_burnin = 100, q = 0.5, params = model_params,
#'                  kernel = "ts")
#'
#' plot(out)
#'
#' @rdname plot.DetectCpObj
#' @export
#'
plot.DetectCpObj <- function(x, y = NULL,
                             plot_freq = FALSE,
                             loss = "VI",
                             maxNClusters = 0,
                             nRuns = 16,
                             maxZealousAttempts = 100, ...) {


  time <- V2 <- y <- obs <- NULL

  if(!plot_freq){

    if(x$kernel_ts){

      if(x$univariate_ts){

        est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                    nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

        cp <- cumsum(table(est_cp))[-length(table(est_cp))] + 1

        vec_data <- x$data

        .data_plot <- as.data.frame(cbind(vec_data))
        .data_plot$time <- rep(1:length(x$data))
        .data_plot$obs <- as.factor(rep(1, nrow(.data_plot)))

        p1 <- ggplot2::ggplot(.data_plot) +
          ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data, color = obs),  linetype = 1) +
          ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
          ggplot2::labs(x = "Time",
                        y = "Value",
                        color = NULL) +
          ggplot2::scale_color_viridis_d() +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position="none")

        print(p1)


      } else {

        est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                    nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

        cp <- cumsum(table(est_cp))[-length(table(est_cp))] + 1

        vec_data <- as.numeric()

        for(i in 1:nrow(x$data)){
          vec_data <- c(vec_data,x$data[i,])
        }

        .data_plot <- as.data.frame(cbind(vec_data, sort(rep(1:nrow(x$data),ncol(x$data)))))
        .data_plot$V2 <- factor(.data_plot$V2, labels = unique(paste0("d = ", .data_plot$V2)) )
        .data_plot$time <- rep(1:ncol(x$data),nrow(x$data))

        p1 <- ggplot2::ggplot(.data_plot) +
          ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data, color = V2),  linetype = 1) +
          ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
          ggplot2::labs(x = "Time",
                        y = "Value",
                        color = NULL) +
          ggplot2::scale_color_viridis_d() +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position="top", legend.key.width = ggplot2::unit(1, 'cm'))

        print(p1)

      }

    }

    if(x$kernel_epi){

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      x$data = t(x$data)

      .df_sf_plot <- data.frame(as.vector(sapply(1:nrow(x$data), function(y) 1 - cumsum(x$data[y,]) / sum(x$data[y,]))),
                                rep(1:ncol(x$data), nrow(x$data)),
                                rep(1:ncol(x$data),nrow(x$data)),
                                rep("1", nrow(x$data)))

      colnames(.df_sf_plot) <- c("y","x","time","obs")

      cp <- cumsum(table(est_cp))[-length(table(est_cp))] + 1


      p1 <- ggplot2::ggplot(.df_sf_plot, ggplot2::aes(x = x, y = y, color = obs)) +
        ggplot2::geom_line(lwd = 0.5) +
        ggplot2::geom_vline(xintercept = unique(.df_sf_plot$time)[cp], linetype = 2) +
        ggplot2::labs(x = "Time",
                      y = "Proportion of Infected Individuals",
                      color = NULL) +
        ggplot2::scale_color_viridis_d() +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="none")

      print(p1)


    }

  } else {

    if(x$kernel_ts){

      if(x$univariate_ts){

        est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                    nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

        cp <- cumsum(table(est_cp))[-length(table(est_cp))] + 1

        vec_data <- x$data

        .data_plot <- as.data.frame(cbind(vec_data))
        .data_plot$time <- 1:length(vec_data)
        .data_plot$obs <- as.factor(rep(1, nrow(.data_plot)))

        p1 <- ggplot2::ggplot(.data_plot) +
          ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data, color = obs),  linetype = 1) +
          ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
          ggplot2::labs(x = " ",
                        y = "Value",
                        color = NULL) +
          ggplot2::scale_colour_brewer(palette = "Set1") +
          ggplot2::theme_minimal()

        p1 <- p1 + ggplot2::theme(legend.position="none")

        x_unique <- unique(.data_plot$time)
        b <- rep(0, length(x$data))

        for(i in 1:x$n_iterations){

          cp_iteration <- cumsum(table(x$orders[i,]))[-length(table(x$orders[i,]))] + 1

          b[cp_iteration] = b[cp_iteration] + 1

        }

        b <- b/(x$n_iterations)

        p2 <- ggplot2::ggplot(data.frame(x = x_unique, y =b)) +
          ggplot2::geom_bar(ggplot2::aes(x = x_unique, y = y), stat="identity", width = 0.5, col = "black") +
          ggplot2::theme_linedraw() +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(angle = 90)) +
          ggplot2::scale_y_continuous(breaks = c(0,.5,1)) +
          ggplot2::ylab("Prob.") +
          ggplot2::xlab("Time") +
          ggplot2::scale_color_viridis_d() +
          ggplot2::theme_minimal()

        p2 <- p2 + ggplot2::theme(legend.position="none")

        print(ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(2,1), common.legend = FALSE))

      } else {

        est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                    nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

        cp <- cumsum(table(est_cp))[-length(table(est_cp))] + 1

        vec_data <- as.numeric()

        for(i in 1:nrow(x$data)){
          vec_data <- c(vec_data,x$data[i,])
        }

        .data_plot <- as.data.frame(cbind(vec_data, sort(rep(1:nrow(x$data),ncol(x$data)))))
        .data_plot$V2 <- factor(.data_plot$V2, labels = unique(paste0("d = ", .data_plot$V2)) )
        .data_plot$time <- rep(1:ncol(x$data),nrow(x$data))

        p1 <- ggplot2::ggplot(.data_plot) +
          ggplot2::geom_line(ggplot2::aes(x = time, y = vec_data, color = V2),  linetype = 1) +
          ggplot2::geom_vline(xintercept = unique(.data_plot$time)[cp], linetype = 2) +
          ggplot2::labs(x = " ",
                        y = "Value",
                        color = NULL) +
          ggplot2::scale_color_viridis_d() +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position="top", legend.key.width = ggplot2::unit(1, 'cm'))


        x_unique <- unique(.data_plot$time)
        b <- rep(0, ncol(x$data))

        for(i in 1:x$n_iterations){

          cp_iteration <- cumsum(table(x$orders[i,]))[-length(table(x$orders[i,]))] + 1

          b[cp_iteration] = b[cp_iteration] + 1

        }

        b <- b/(x$n_iterations)

        p2 <- ggplot2::ggplot(data.frame(x = x_unique, y =b)) +
          ggplot2::geom_bar(ggplot2::aes(x = x_unique, y = y), stat="identity", width = 0.5, col = "black") +
          ggplot2::theme_linedraw() +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(angle = 90)) +
          ggplot2::scale_y_continuous(breaks = c(0,.5,1)) +
          ggplot2::ylab("Prob.") +
          ggplot2::xlab("Time") +
          ggplot2::theme_minimal()

        print(ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(2,1), common.legend = TRUE))

      }

    }

    if(x$kernel_epi){

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      x$data = t(x$data)

      .df_sf_plot <- data.frame(as.vector(sapply(1:nrow(x$data), function(y) 1 - cumsum(x$data[y,]) / sum(x$data[y,]))),
                                rep(1:ncol(x$data), nrow(x$data)),
                                rep(1:ncol(x$data),nrow(x$data)),
                                rep("1", nrow(x$data)))

      colnames(.df_sf_plot) <- c("y","x","time","obs")

      cp <- cumsum(table(est_cp))[-length(table(est_cp))] + 1


      p1 <- ggplot2::ggplot(.df_sf_plot, ggplot2::aes(x = x, y = y, color = obs)) +
        ggplot2::geom_line(lwd = 0.5) +
        ggplot2::geom_vline(xintercept = unique(.df_sf_plot$time)[cp], linetype = 2) +
        ggplot2::labs(x = "Time",
                      y = "Proportion of Infected Individuals",
                      color = NULL) +
        ggplot2::scale_color_viridis_d() +
        ggplot2::theme_minimal()

      p1 <- p1 + ggplot2::theme(legend.position="none")

      x_unique <- unique(.df_sf_plot$time)
      b <- rep(0, length(x$data))

      for(i in 1:x$n_iterations){

        cp_iteration <- cumsum(table(x$orders[i,]))[-length(table(x$orders[i,]))] + 1

        b[cp_iteration] = b[cp_iteration] + 1

      }

      b <- b/(x$n_iterations)

      p2 <- ggplot2::ggplot(data.frame(x = x_unique, y =b)) +
        ggplot2::geom_bar(ggplot2::aes(x = x_unique, y = y), stat="identity", width = 0.5, col = "black") +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(angle = 90)) +
        ggplot2::scale_y_continuous(breaks = c(0,.5,1)) +
        ggplot2::ylab("Prob.") +
        ggplot2::xlab("Time") +
        ggplot2::scale_color_viridis_d() +
        ggplot2::theme_minimal()

      p2 <- p2 + ggplot2::theme(legend.position="none")

      print(ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(2,1), common.legend = FALSE))

    }

  }

}
