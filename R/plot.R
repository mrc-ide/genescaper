
#------------------------------------------------
#' @title Produce visualisation of prior on allele frequencies
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param n_reps TODO
#' @param p_step TODO
#'
#' @importFrom stats rgamma rnorm
#' @export

plot_prior_freq <- function(project, n_reps = 10, p_step = 0.01) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(n_reps, zero_allowed = FALSE)
  assert_leq(n_reps, 1e3)
  assert_single_bounded(p_step, right = 0.1, inclusive_left = FALSE)
  
  # get model parameters
  mu_mean <- project$model$parameters$mu_mean
  mu_scale <- project$model$parameters$mu_scale
  sigsq_mean <- project$model$parameters$sigsq_mean
  sigsq_var <- project$model$parameters$sigsq_var
  
  # reparameterise for convenience
  sigsq_shape <- sigsq_mean^2 / sigsq_var + 2
  sigsq_scale <- sigsq_mean*(sigsq_shape - 1)
  
  # draw from normal-inverse-gamma prior (note that rate of gamma becomes scale
  # of inverse gamma)
  sigsq <- 1.0 / rgamma(n_reps, shape = sigsq_shape, rate = sigsq_scale)
  mu <- rnorm(n_reps, mean = mu_mean, sd = sqrt(sigsq * mu_scale))
  
  # get logit-normal distribution over all reps
  p <- seq(p_step, 1.0 - p_step, p_step)
  fp <- mapply(function(i) {
    dlogitnorm(p, raw_mu = mu[i], raw_sigsq = sigsq[i], return_log = FALSE)
  }, 1:n_reps)
  
  # create plot object
  plot1 <- expand.grid( p = p, rep = 1:n_reps) %>%
    dplyr::bind_cols(fp = as.vector(fp)) %>%
    ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(x = p, y = fp, group = rep)) +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::xlab("Allele frequency") + ggplot2::ylab("Probability")
  
  return(plot1)
}

#------------------------------------------------
#' @title Produce visualisation of prior on spatial correlation
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param type TODO
#' @param n_reps TODO
#' @param d_max TODO
#'
#' @importFrom stats rgamma rbeta
#' @export

plot_prior_cor <- function(project, type = 1, n_reps = 10, d_max = NULL) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(n_reps, zero_allowed = FALSE)
  assert_leq(n_reps, 1e3)
  assert_in(type, 1:2)
  if (!is.null(d_max)) {
    assert_single_pos(d_max, zero_allowed = FALSE)
  }
  
  # get model parameters
  lambda_mean <- project$model$parameters$lambda_mean
  lambda_var <- project$model$parameters$lambda_var
  nu_shape1 <- project$model$parameters$nu_shape1
  nu_shape2 <- project$model$parameters$nu_shape2
  
  # reparameterise for convenience
  lambda_shape <- lambda_mean^2 / lambda_var
  lambda_rate <- lambda_mean / lambda_var
  
  # split on type
  if (type == 1) {
    # plot ribbon
    
  } else if (type == 2) {
    # plot random draws
    
    # draw from priors
    lambda <- rgamma(n_reps, shape = lambda_shape, rate = lambda_rate)
    nu <- rbeta(n_reps, shape1 = nu_shape1, shape2 = nu_shape2)
    
    # choose d_max automatically from prior parameters
    if (is.null(d_max)) {
      d_max <- 3 * max(lambda)
    }
    
    # get correlation curve over all reps
    d <- seq(0, d_max, l = 201)
    fd <- mapply(function(i) {
      (1.0 - nu[i])*exp(-d / lambda[i])
    }, 1:n_reps)
    
    # create plot object
    plot1 <- expand.grid(d = d, rep = 1:n_reps) %>%
      dplyr::bind_cols(fd = as.vector(fd)) %>%
      ggplot2::ggplot() + ggplot2::theme_bw() +
      ggplot2::geom_line(ggplot2::aes(x = d, y = fd, group = rep)) +
      ggplot2::ylim(c(0, 1)) +
      ggplot2::xlab("Distance") + ggplot2::ylab("Correlation")
  }
  
  return(plot1)
}
