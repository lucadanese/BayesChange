#' ClustCpObj class constructor
#'
#' @description
#' A constructor for the \code{ClustCpObj} class, which stores the output of
#' the change point detection and clustering algorithms.
#'
#' @param data A vector or matrix containing the observed data.
#' @param n_iterations Total number of MCMC iterations.
#' @param n_burnin Number of burn-in iterations removed from posterior summaries.
#' @param clust A matrix where each row contains the cluster assignments for one iteration.
#' @param orders A multidimensional array where each slice is a matrix representing
#'   the latent order at each iteration.
#' @param time Total computational time (in seconds).
#' @param norm_vec A vector containing precomputed normalization constants.
#' @param entropy_MCMC A \code{coda::mcmc} object containing the MCMC samples of the entropy.
#' @param lkl_MCMC A \code{coda::mcmc} object containing the log-likelihood values at each iteration.
#' @param I0_MCMC A \code{coda::mcmc} object with the MCMC trace of the initial infection proportion \eqn{I_0}.
#' @param kernel_ts Logical; TRUE if the kernel corresponds to time-series data.
#' @param kernel_epi Logical; TRUE if the kernel corresponds to epidemic diffusion data.
#' @param univariate_ts Logical; TRUE if the data represent a univariate time series,
#'   FALSE for multivariate time series.
#'
#' @return An object of class \code{ClustCpObj}.
#' @export
#'
ClustCpObj <- function(data = NULL,
                       n_iterations = NULL,
                       n_burnin = NULL,
                       clust = NULL,
                       orders = NULL,
                       time = NULL,
                       norm_vec = NULL,
                       entropy_MCMC = NULL,
                       lkl_MCMC = NULL,
                       I0_MCMC = NULL,
                       kernel_ts = NULL,
                       kernel_epi = NULL,
                       univariate_ts = NULL){

  value <- list(data = data,
                n_iterations = n_iterations,
                n_burnin = n_burnin,
                clust = clust,
                orders = orders,
                time = time,
                norm_vec = norm_vec,
                entropy_MCMC = entropy_MCMC,
                lkl_MCMC = lkl_MCMC,
                I0_MCMC = I0_MCMC,
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
#' data("stock_uni")
#'
#' params_uni <- list(a = 1,
#'                    b = 1,
#'                    c = 1,
#'                    phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:3,], n_iterations = 1000, n_burnin = 100,
#'                 L = 1, q = 0.5, B = 500, params = params_uni, kernel = "ts")
#'
#' print(out)
#'
#' @rdname print.ClustCpObj
#' @export
#'
print.ClustCpObj <- function(x, ...) {
  cat("ClustCpObj object\n")
  if(x$kernel_ts){
    if(x$univariate_ts){
      cat("Type: clustering univariate time series with common change points\n")
    } else {
      cat("Type: clustering multivariate time series with common change points\n")
    }
  } else if(x$kernel_epi){
    cat("Type: clustering epidemic diffusions with common change points\n")
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
#' data("stock_uni")
#'
#' params_uni <- list(a = 1,
#'                    b = 1,
#'                    c = 1,
#'                    phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:3,], n_iterations = 1000, n_burnin = 100,
#'                 L = 1, q = 0.5, B = 500, params = params_uni, kernel = "ts")
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
      cat("Clustering summary (univariate time series):\n",
          "- Burn-in iterations:", object$n_burnin, "\n",
          "- MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
          "- Average number of clusters:",
          round(mean(apply(object$clust[(object$n_burnin + 1):object$n_iterations, ], 1, function(x) length(unique(x)))), 2), "\n",
          "- Computational time:", round(object$time, 2), "seconds\n",
          "\nUse plot() for a detailed visualization or posterior_estimate() to analyze the clustering results.\n")
    } else {
      cat("Clustering ", paste0(nrow(object$data[,,1]), "-dimensional time series:\n"),
          "- Burn-in iterations:", object$n_burnin, "\n",
          "- MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
          "- Average number of clusters:",
          round(mean(apply(object$clust[(object$n_burnin + 1):object$n_iterations, ], 1, function(x) length(unique(x)))), 2), "\n",
          "- Computational time:", round(object$time, 2), "seconds\n",
          "\nUse plot() for a detailed visualization or posterior_estimate() to analyze the clustering results.\n")
    }
  } else if (object$kernel_epi){
    cat("Clustering epidemic diffusions:\n",
        "- Burn-in iterations:", object$n_burnin, "\n",
        "- MCMC iterations:", object$n_iterations - object$n_burnin, "\n",
        "- Average number of clusters:",
        round(mean(apply(object$clust[(object$n_burnin + 1):object$n_iterations, ], 1, function(x) length(unique(x)))), 2), "\n",
        "- Computational time:", round(object$time, 2), "seconds\n",
        "\nUse plot() for a detailed visualization or posterior_estimate() to analyze the clustering results.\n")
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
#' data("stock_uni")
#'
#' params_uni <- list(a = 1,
#'                    b = 1,
#'                    c = 1,
#'                    phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:3,], n_iterations = 1000, n_burnin = 100,
#'                 L = 1, q = 0.5, B = 500, params = params_uni, kernel = "ts")
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
#' data("stock_uni")
#'
#' params_uni <- list(a = 1,
#'                    b = 1,
#'                    c = 1,
#'                    phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:3,], n_iterations = 1000, n_burnin = 100,
#'                 L = 1, q = 0.5, B = 500, params = params_uni, kernel = "ts")
#'
#' plot(out)
#'
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
      #.data_plot <- dplyr::mutate(.data_plot, Observation = as.numeric(gsub("V", "", Observation)))
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
        ggplot2::scale_color_viridis_d() +
        ggplot2::theme_minimal()

    } else {

      est_cp = posterior_estimate(x, loss = loss, maxNClusters = maxNClusters,
                                  nRuns = nRuns, maxZealousAttempts = maxZealousAttempts)

      .data_plot <- data.frame(Value = numeric(0))

      obs_names <- dimnames(x$data)[[3]]

      count = 0
      for (i in 1:dim(x$data)[3]) {
        mat <- x$data[,,i]

        obs_label <- if (!is.null(obs_names)) obs_names[i] else i

        for (j in 1:nrow(mat)) {
          count = count + 1

          .data_plot <- rbind(
            .data_plot,
            data.frame(
              Value = as.vector(mat[j, ]),
              Observation = obs_label,
              Cluster = est_cp[i],
              Time = 1:ncol(mat),
              Count = count
            )
          )
        }
      }

      .data_plot$Observation <- as.factor(.data_plot$Observation)
      .data_plot$Cluster <- as.factor(.data_plot$Cluster)
      .data_plot$Count <- as.factor(.data_plot$Count)

      ggplot2::ggplot(.data_plot) +
        ggplot2::geom_line(
          ggplot2::aes(
            x = Time,
            y = Value,
            color = Observation,
            group = Count,
            linetype = Cluster
          )
        ) +
        ggplot2::xlab("Time") +
        ggplot2::ylab("Value") +
        ggplot2::scale_color_viridis_d() +
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
      ggplot2::scale_color_viridis_d() +
      ggplot2::theme_minimal()

  }

}

#' set generic
#' @name posterior_estimate
#' @keywords internal
#' @export
#'
plot_psm <- function (object, ...) {
  UseMethod("plot_psm")
}

#' Plot the Posterior Similarity Matrix (PSM) for a ClustCpObj
#'
#' This function computes and visualizes the posterior similarity matrix
#' (PSM) from a \code{ClustCpObj} object.
#' The PSM shows the posterior co-clustering probabilities of all observations.
#'
#' @param object an object of class \code{ClustCpObj}.
#' @param reorder Logical; if \code{TRUE} (default), items are reordered using
#'   hierarchical clustering to highlight clusters in the final plot
#' @param title Character; the plot title (default: \code{"Posterior Similarity Matrix"}).
#' @param ... parameter of the generic method.
#'
#' @return A \code{ggplot2} object representing the posterior similarity matrix.
#'
#' @importFrom rlang .data
#'
#' @rdname plot_psm.ClustCpObj
#' @examples
#' data("stock_uni")
#'
#' params_uni <- list(a = 1,
#'                    b = 1,
#'                    c = 1,
#'                    phi = 0.1)
#'
#' out <- clust_cp(data = stock_uni[1:3,], n_iterations = 1000, n_burnin = 100,
#'                 L = 1, q = 0.5, B = 500, params = params_uni, kernel = "ts")
#' plot_psm(out)
#' @export
#'
plot_psm.ClustCpObj <- function(object,
                                reorder = TRUE,
                                title = "Posterior Similarity Matrix", ...) {
  if (!inherits(object, "ClustCpObj")) {
    stop("object must be of class 'ClustCpObj'")
  }

  # Extract MCMC cluster samples
  if (is.null(object$clust)) {
    stop("The ClustCpObj does not contain MCMC cluster samples (clust).")
  }

  mcmc_chain <- object$clust[(object$n_burnin + 1):object$n_iterations,]
  n_items <- ncol(mcmc_chain)

  # Compute Posterior Similarity Matrix using SALSO
  psm <- salso::psm(mcmc_chain)

  obs_labels <- 1:n_items

  # Reorder items if requested
  if (reorder) {
    hc <- stats::hclust(stats::dist(1 - psm))
    order_idx <- hc$order
    psm <- psm[order_idx, order_idx]
    obs_labels <- obs_labels[order_idx]
  }

  psm_melt <- reshape2::melt(psm[1:nrow(psm), , drop = FALSE])
  colnames(psm_melt) <- c("Obs1", "Obs2", "Similarity")

  # Plot
  plot_heat <- ggplot2::ggplot(psm_melt) +
    ggplot2::geom_tile(
      ggplot2::aes(x = .data$Obs1, y = .data$Obs2, fill = .data$Similarity),
      na.rm = TRUE,
      linewidth = 0
    ) +
    ggplot2::scale_fill_gradient(low = "transparent", high = "#470D60FF") +
    ggplot2::labs(
      title = title,
      fill = "Co-clustering\nProbability",
      x = " ",
      y = " "
    ) +
    ggplot2::theme_linedraw() +
    ggplot2::coord_cartesian(
      xlim = c(0, ncol(psm) + 1),
      ylim = c(0, ncol(psm) + 1),
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, ncol(psm), by = 1),
      labels = obs_labels
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(1, ncol(psm), by = 1),
      labels = obs_labels
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_line(linewidth = 0.1)
    )

  return(plot_heat)
}
