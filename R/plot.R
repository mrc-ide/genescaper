
#------------------------------------------------
#' @title Red to blue colours
#'
#' @description Simple sequence of red-to-blue colours.
#'
#' @param n the number of colours.
#'
#' @importFrom grDevices colorRampPalette
#' @export

col_hotcold <- function(n = 6) {
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  my_pal <- colorRampPalette(raw_cols)
  return(my_pal(n))
}

#------------------------------------------------
#' @title Plot pairwise spatial distance vs. genetic distances
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param genetic_stat TODO
#'
#' @export

plot_pairwise <- function(project, genetic_stat = "Gst") {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_in(genetic_stat, c("Gst", "D"))
  
  # get spatial and genetic distance matrices
  spatial_dist <- project$data$pairwise_measures$distance
  if (genetic_stat == "Gst") {
    genetic_dist <- project$data$pairwise_measures$Gst
  } else if (genetic_stat == "D") {
    genetic_dist <- project$data$pairwise_measures$D
  }
  
  # get into dataframe for plotting
  df_plot <- data.frame(spatial = spatial_dist[upper.tri(spatial_dist, diag = FALSE)],
                        genetic = genetic_dist[upper.tri(genetic_dist, diag = FALSE)])
  
  # produce plot
  ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(x = .data$spatial, y = .data$genetic)) +
    ggplot2::xlab("spatial distance (km)") +
    ggplot2::ylab(sprintf("pairwise %s", genetic_stat))
  
}
#------------------------------------------------
#' @title Plot prior and/or posterior spatial correlation
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param prior_on TODO
#' @param posterior_on TODO
#' @param n_reps TODO
#' @param d_max TODO
#'
#' @importFrom stats rgamma rbeta runif
#' @importFrom grDevices grey
#' @export

plot_spatial_autocor <- function(project, prior_on = TRUE, posterior_on = TRUE, n_reps = 1e3, d_max = NULL) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_logical(prior_on)
  assert_single_logical(posterior_on)
  assert_single_pos_int(n_reps, zero_allowed = FALSE)
  assert_leq(n_reps, 1e6)
  if (!is.null(d_max)) {
    assert_single_pos(d_max, zero_allowed = FALSE)
  }
  if (!prior_on && !posterior_on) {
    stop("at least one of prior_on and posterior_on must be TRUE")
  }
  
  # get prior draws
  if (prior_on) {
    model_params <- project$model$parameters
    prior_draws <- data.frame(lambda = rgamma(n_reps, shape = model_params$lambda_shape, rate = model_params$lambda_rate),
                              nu = rbeta(n_reps, shape1 = model_params$nu_shape1, shape2 = model_params$nu_shape2),
                              omega = runif(n_reps, min = 1, max = 3))
  }
  
  # establish whether posterior draws exist in project
  posterior_exists <- !is.null(project$model$MCMC$output)
  
  # get posterior draws
  if (posterior_exists) {
    post_draws <- project$model$MCMC$output %>%
      dplyr::filter(.data$phase == "sampling") %>%
      dplyr::select(.data$lambda, .data$nu, .data$omega) %>%
      dplyr::sample_n(n_reps, replace = TRUE)
  }
  
  # if plotting posterior then choose default d_max based on posterior
  if (is.null(d_max) && posterior_on && posterior_exists) {
    d_max <- post_draws %>%
      dplyr::mutate(max_dist = .data$lambda * (-log(0.01))^(1 / .data$omega)) %>%
      dplyr::pull(.data$max_dist) %>%
      max()
  }
  
  # if not plotting posterior then choose default d_max based on prior
  if (is.null(d_max)) {
    d_max <- prior_draws %>%
      dplyr::mutate(max_dist = .data$lambda * (-log(0.01))^(1 / .data$omega)) %>%
      dplyr::pull(.data$max_dist) %>%
      max()
  }
  
  # define plotting distances
  d <- seq(0, d_max, l = 101)
  
  # define helper function for creating function with distance
  f_d <- function(p) {
    p[2] * exp(-(d / p[1])^p[3])
  }
  
  # generate prior bubble
  if (prior_on) {
    prior_mat <- apply(prior_draws, 1, f_d)
    prior_quants <- apply(prior_mat, 1, function(x) quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))) %>%
      t() %>% as.data.frame()
    names(prior_quants) <- sprintf("Q%s", c(2.5, 10, 50, 90, 97.5))
  }
  
  # generate posterior bubble
  if (posterior_on && posterior_exists) {
    post_mat <- apply(post_draws, 1, f_d)
    post_quants <- apply(post_mat, 1, function(x) quantile(x, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))) %>%
      t() %>% as.data.frame()
    names(post_quants) <- sprintf("Q%s", c(2.5, 10, 50, 90, 97.5))
  }
  
  # create empty plot
  plot1 <- ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::xlim(c(0, d_max)) + ggplot2::ylim(c(0,1)) +
    ggplot2::xlab("spatial distance (km)") +
    ggplot2::ylab("autocorrelation")
  
  # add prior bubble
  if (prior_on) {
    plot1 <- plot1 +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Q2.5, ymax = .data$Q97.5, x = d),
                           fill = grey(0.5), alpha = 0.2, data = prior_quants) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Q10, ymax = .data$Q90, x = d),
                           fill = grey(0.5), alpha = 0.2, data = prior_quants) +
      ggplot2::geom_line(ggplot2::aes(x = d, y = .data$Q50), col = grey(0.5), data = prior_quants)
  }
  
  # add posterior bubble
  if (posterior_on && posterior_exists) {
    plot1 <- plot1 +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Q2.5, ymax = .data$Q97.5, x = d),
                           fill = "red", alpha = 0.2, data = post_quants) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Q10, ymax = .data$Q90, x = d),
                           fill = "red", alpha = 0.2, data = post_quants) +
      ggplot2::geom_line(ggplot2::aes(x = d, y = .data$Q50), col = "red", data = post_quants)
  }
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot GeoMAPI coverage
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param breaks the sequence of coverage breaks used.
#' @param col_scale the colour scale to use.
#' @param plot_sampling_points whether to overlay sampling locations.
#' 
#' @import sf
#' @import ggplot2
#' @export

GeoMAPI_plot_coverage <- function(project, breaks = c(0, 10, 20, 30, 40, 50, 100, Inf),
                                  col_scale = col_hotcold, plot_sampling_points = TRUE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_vector_pos(breaks)
  assert_increasing(breaks)
  assert_greq(length(breaks), 2)
  assert_class(col_scale, "function")
  assert_single_logical(plot_sampling_points)
  
  # bin coverage
  coverage <- project$GeoMAPI$coverage
  intersect_bin <- cut(coverage, breaks = breaks, right = FALSE)
  
  # basic plot
  plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank())
  
  # add polys
  plot1 <- plot1 + geom_sf(aes_(fill = ~intersect_bin), color = NA, data = project$maps$grid$polygons)
  
  # add points
  if (plot_sampling_points) {
    plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude),
                                data = project$data$raw$site_data, size = 0.5)
  }
  
  # titles and legends
  plot1 <- plot1 + scale_fill_manual(values = col_scale(length(breaks) - 1),
                                     limits = levels(intersect_bin),
                                     name = "Coverage")
  plot1 <- plot1 + xlab("Longitude") + ylab("Latitude")
  plot1 <- plot1 + guides(fill = guide_legend(reverse = TRUE))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot GeoMAPI z-scores and significance
#'
#' @description TODO
#'
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param plot_sampling_points whether to overlay sampling locations.
#' @param plot_significance whether to outline areas that were identified as
#'   significant outliers.
#' @param min_hex_coverage minimum coverage (number of edges assigned to a hex)
#'   for it to be included in the final result, otherwise these hexes are given
#'   the value \code{NA}.
#' @param col_scale the colour scale to use.
#' @param zlim the limits of the colour scale. If \code{NULL} then these limits
#'   are chosen automatically.
#' @param base_plot optional base plot (object of class \code{ggplot}) on which
#'   this function builds. If \code{NULL} then a simple empty plot is used.
#' @param point_size,point_colour,point_fill,point_stroke properties of plotted
#'   sampling points.
#' 
#' @import sf
#' @import ggplot2
#' @export

GeoMAPI_plot_zscore <- function(project,
                                plot_sampling_points = TRUE,
                                plot_significance = TRUE,
                                min_hex_coverage = 10,
                                col_scale = viridisLite::viridis,
                                zlim = NULL,
                                base_plot = NULL,
                                point_size = 1,
                                point_colour = "white",
                                point_fill = "black",
                                point_stroke = 0.2) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_logical(plot_sampling_points)
  assert_single_logical(plot_significance)
  assert_single_pos_int(min_hex_coverage, zero_allowed = TRUE)
  assert_class(col_scale, "function")
  if (!is.null(zlim)) {
    assert_limit(zlim)
  }
  if (!is.null(base_plot)) {
    assert_class(base_plot, "ggplot")
  }
  assert_single_pos(point_size, zero_allowed = FALSE)
  assert_single_string(point_colour)
  assert_single_string(point_fill)
  assert_single_pos(point_stroke, zero_allowed = FALSE)
  
  # get z_scores and replace with NA if below minimum coverage 
  z_score <- project$GeoMAPI$z_score
  z_score[project$GeoMAPI$coverage < min_hex_coverage] <- NA
  
  # produce basic plot
  if (is.null(base_plot)) {
    plot1 <- ggplot() + theme_bw() + theme(panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank())
  } else {
    plot1 <- base_plot
  }
  
  # add z-scores
  plot1 <- plot1 + geom_sf(aes_(fill = ~z_score), color = NA, data = project$maps$grid$polygons)
  
  # outline significance
  if (plot_significance) {
    
    # outline high values
    w_upper <- project$GeoMAPI$significance$upper
    if (length(w_upper) != 0) {
      merged_poly_upper <- get_merged_poly(project$maps$grid$polygons[w_upper],
                                           d = project$maps$grid$parameters$hex_width / 10)
      plot1 <- plot1 + geom_sf(color = "white", fill = NA, data = merged_poly_upper)
    }
    
    # outline low values
    w_lower <- project$GeoMAPI$significance$lower
    if (length(w_lower) != 0) {
      merged_poly_lower <- get_merged_poly(project$maps$grid$polygons[w_lower],
                                           d = project$maps$grid$parameters$hex_width / 10)
      plot1 <- plot1 + geom_sf(color = "black", fill = NA, data = merged_poly_lower)
    }
    
  }
  
  
  # add sampline points
  if (plot_sampling_points) {
    plot1 <- plot1 + geom_point(aes_(x = ~longitude, y = ~latitude),
                                shape = 21, color = point_colour, fill = point_fill, size = point_size,
                                stroke = point_stroke, data = project$data$raw$site_data)
  }
  
  # titles and legends
  plot1 <- plot1 + scale_fill_gradientn(colours = col_scale(100), name = "z-score", limits = zlim)
  plot1 <- plot1 + xlab("Longitude") + ylab("Latitude")
  
  # return plot object
  return(plot1)
}
