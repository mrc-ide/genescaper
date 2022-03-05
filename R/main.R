#------------------------------------------------
#' @title Load data into project
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param site_data TODO
#' @param genetic_data TODO
#'
#' @importFrom rlang .data
#' @export

bind_data <- function(project, site_data, genetic_data) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  
  # check format of site_data
  assert_dataframe(site_data)
  assert_in(c("site_ID", "longitude", "latitude"), names(site_data),
            message = "site_data must have columns {site_ID, longitude, latitude}")
  assert_vector_numeric(site_data$longitude)
  assert_bounded(site_data$longitude, left = -180, right = 180)
  assert_vector_numeric(site_data$latitude)
  assert_bounded(site_data$latitude, left = -90, right = 90)
  
  # check format of genetic_data
  assert_dataframe(genetic_data)
  assert_in(c("site_ID", "locus", "allele", "freq"), names(genetic_data),
            message = "genetic_data must have columns {site_ID, locus, allele, freq}")
  assert_vector_pos_int(genetic_data$locus)
  assert_vector_pos_int(genetic_data$allele)
  assert_vector_pos(genetic_data$freq)
  assert_bounded(genetic_data$freq, inclusive_left = FALSE, inclusive_right = FALSE,
                 message = "allele frequencies must be in the range (0, 1), and cannot equal exactly 0 or 1")
  
  # check that site_IDs match between datasets
  assert_eq(sort(unique(site_data$site_ID)),
            sort(unique(genetic_data$site_ID)),
            message = "site_data and genetic_data must contain the same set of site_ID values")
  
  # check that all loci are represented in all sites and alleles
  locus_match <- genetic_data %>%
    dplyr::group_by(.data$site_ID, .data$allele) %>%
    dplyr::summarise(hash = rlang::hash(.data$locus),
                     .groups = "drop_last")
  if (length(unique(locus_match$hash)) != 1) {
    stop("the same set of loci (in the same order) must be represented for all site_ID & allele combinations.")
  }
  
  # check that the same alleles are represented in all sites and loci
  allele_match <- genetic_data %>%
    dplyr::group_by(.data$site_ID, .data$locus) %>%
    dplyr::summarise(hash = rlang::hash(.data$locus),
                     .groups = "drop_last") %>%
    dplyr::group_by(.data$locus) %>%
    dplyr::summarise(same_alleles = length(unique(.data$hash)) == 1)
  if (!all(allele_match$same_alleles)) {
    stop("for a given locus, the same set of alleles (in the same order) must be represented for all site_IDs.")
  }
  
  # load data into project
  project$data$raw <- list(site_data = site_data,
                           genetic_data = genetic_data)
  
  # get pairwise spatial (great circle) distance between sites and load into
  # project
  project$data$pairwise_measures$distance <- get_GC_distance(site_data$longitude, site_data$latitude)
  
  return(project)
}

#------------------------------------------------
#' @title Calculate pairwise Gst between sites
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#'
#' @export

get_pairwise_Gst <- function(project) {
  
  freq_array <- project$data$raw$genetic_data %>%
    dplyr::group_by(.data$locus) %>%
    dplyr::group_split() %>%
    lapply(function (x) {
      x %>% dplyr::group_by(.data$allele) %>%
        dplyr::group_split() %>%
        lapply(function (x) x$freq)
    })
  
  # get allele frequencies into list over loci then matrix over demes and
  # alleles
  freq_list <- project$data$raw$genetic_data %>%
    dplyr::group_by(.data$locus) %>%
    dplyr::group_split() %>%
    lapply(function (x) {
      x %>% tidyr::pivot_wider(names_from = .data$allele, values_from = .data$freq) %>%
        dplyr::select(-.data$site_ID, -.data$locus) %>%
        as.matrix()
    })
  
  # pass to efficient C++ function
  output_raw <- get_mean_pairwise_Gst_cpp(freq_list)
  
  # convert to distance matrix and load into project
  project$data$pairwise_measures$Gst <- as.dist(t(output_raw))
  
  return(project)
}

#------------------------------------------------
#' @title Define parameters of the spatial correlation model
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param mu_mean TODO
#' @param mu_scale TODO
#' @param sigsq_mean TODO
#' @param sigsq_var TODO
#' @param lambda_mean TODO
#' @param lambda_var TODO
#' @param nu_shape1 TODO
#' @param nu_shape2 TODO
#'
#' @export

define_model <- function(project, mu_mean = 0.0, mu_scale = 1.0, sigsq_mean = 1.0, sigsq_var = 1.0,
                         lambda_mean = 5.0, lambda_var = 10.0, nu_shape1 = 1.0, nu_shape2 = 1.0) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_numeric(mu_mean)
  assert_single_pos(mu_scale, zero_allowed = FALSE)
  assert_single_pos(sigsq_mean, zero_allowed = FALSE)
  assert_single_pos(sigsq_var, zero_allowed = FALSE)
  assert_single_pos(lambda_mean, zero_allowed = FALSE)
  assert_single_pos(lambda_var, zero_allowed = FALSE)
  assert_single_pos(nu_shape1, zero_allowed = FALSE)
  assert_single_pos(nu_shape2, zero_allowed = FALSE)
  
  # store within project
  project$model$parameters <- list(mu_mean = mu_mean,
                                   mu_scale = mu_scale,
                                   sigsq_mean = sigsq_mean,
                                   sigsq_var = sigsq_var,
                                   lambda_mean = lambda_mean,
                                   lambda_var = lambda_var,
                                   nu_shape1 = nu_shape1,
                                   nu_shape2 = nu_shape2)
  
  return(project)
}

#------------------------------------------------
#' @title Run MCMC
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param ... additional parameters that will be passed to
#'   \code{drjacoby::run_mcmc()}.
#'
#' @export

run_mcmc <- function(project, ...) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  
  # check that both data and model have been defined
  assert_non_null(project$dat$raw)
  assert_non_null(project$model$parameters)
  
  # get distance between all sites
  site_dist <- as.matrix(project$data$pairwise_measures$distance)
  
  # Transform allele frequencies to z values. These are split into a list over
  # all locus&allele combos. Each combo has an independent mean and variance
  # under the model and so we can treat these as equivalent replicates.
  z_list <- project$data$raw$genetic_data %>%
    transform_p_to_z() %>%
    dplyr::select(.data$locus, .data$allele, .data$z) %>%
    dplyr::group_by(.data$locus, .data$allele) %>%
    dplyr::group_split() %>%
    lapply(function (x) x$z)
  names(z_list) <- seq_along(z_list)
  
  # define parameters dataframe
  df_params <- drjacoby::define_params(name = "log_lambda", min = -Inf, max = Inf,
                                       name = "nu", min = 0, max = 1)
  
  # source C++ likelihood and prior functions
  Rcpp::sourceCpp(system.file("extdata/GRF_model.cpp", package = 'genescaper', mustWork = TRUE))
  
  # run MCMC
  mcmc <- drjacoby::run_mcmc(data = z_list,
                             df_params = df_params,
                             misc = append(project$model$parameters,
                                           list(site_dist = site_dist,
                                                n_site = nrow(site_dist))),
                             loglike = "loglike",
                             logprior = "logprior",
                             ...)
  
  # replace log_lambda with lambda in output
  mcmc$output <- mcmc$output %>%
    dplyr::rename(lambda = .data$log_lambda) %>%
    dplyr::mutate(lambda = exp(.data$lambda))
  
  # add to project
  project$model$MCMC <- mcmc
  
  return(project)
}

#------------------------------------------------
#' @title Create map composed of hexagonal tiles
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param hex_width width of hexagons.
#' @param border_coords dataframe giving coordinates (longitude, latitude) of a
#'   polygon within which the map is defined. If null then this is generated
#'   automatically from the convex hull of the sample site locations.
#'
#' @import sf
#' @importFrom grDevices chull
#' @export

create_hex_grid <- function(project, hex_width = NULL, border_coords = NULL) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  if (!is.null(hex_width)) {
    assert_single_pos(hex_width, zero_allowed = FALSE)
  }
  if (!is.null(border_coords)) {
    assert_dataframe(border_coords)
    assert_in(c("longitude", "latitude"), names(border_coords))
    assert_vector_numeric(border_coords$longitude)
    assert_vector_numeric(border_coords$latitude)
    assert_bounded(border_coords$longitude, left = -180, right = 180)
    assert_bounded(border_coords$latitude, left = -90, right = 90)
  }
  
  message("Creating hex map")
  
  # calculate default hex size from data
  if (is.null(hex_width)) {
    if (is.null(project$data$raw$site_data)) {
      stop("hex_width must be speficied when data is not loaded")
    }
    diff_longitude <- diff(range(project$data$raw$site_data$longitude))
    diff_latitude <- diff(range(project$data$raw$site_data$latitude))
    hex_width <- min(diff_longitude, diff_latitude) / 20
    message(sprintf("hex width chosen automatically: %s", signif(hex_width, 3)))
  }
  
  # get border_coords from convex hull of data
  if (is.null(border_coords)) {
    data_coords <- project$data$raw$site_data %>%
      dplyr::select(.data$longitude, .data$latitude)
    ch_data <- chull(data_coords)
    border_coords <- data_coords[c(ch_data, ch_data[1]),]
  }
  
  # get convex hull into sf polygon format
  bounding_poly <- sf::st_sfc(sf::st_polygon(list(as.matrix(border_coords))))
  
  # make sf hex grid from poly
  hex_polys <- sf::st_make_grid(bounding_poly, cellsize = hex_width, square = FALSE)
  nhex <- length(hex_polys)
  
  # get hex centroid points
  hex_pts <- sf::st_centroid(hex_polys)
  hex_pts_df <- as.data.frame(t(mapply(as.vector, hex_pts)))
  names(hex_pts_df) <- c("longitude", "latitude")
  
  message(sprintf("%s hexagons created", nhex))
  
  # add to project
  project$maps <- list(grid = list(parameters = list(hex_width = hex_width),
                                   centroids = hex_pts_df,
                                   polygons = hex_polys))
  
  return(project)
}

#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param loci TODO
#' @param reps TODO
#' @param inner_reps TODO
#' @param quantiles TODO
#' @param exceedance TODO
#' @param pb_markdown TODO
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats quantile sd
#' @export

predict_map <- function(project, loci, reps = 2, inner_reps = 10,
                        quantiles = c(0.025, 0.5, 0.975), exceedance = NULL,
                        pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_vector_pos_int(loci, zero_allowed = FALSE)
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_single_pos_int(inner_reps, zero_allowed = FALSE)
  if (!is.null(quantiles)) {
    assert_vector_bounded(quantiles)
  }
  if (!is.null(exceedance)) {
    assert_vector_bounded(exceedance)
  }
  assert_single_logical(pb_markdown)
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$lambda, .data$nu) %>%
    dplyr::sample_n(reps) %>%
    as.list()
  
  # transform frequencies to continuous scale
  data_df <- project$data$raw$genetic_data %>%
    dplyr::filter(.data$locus %in% loci) %>%
    transform_p_to_z()
  
  # get z values split by locus and grouped into matrix by allele
  data_list <- lapply(split(data_df, data_df$locus), function(x) {
    tidyr::pivot_wider(x, names_from = .data$site_ID, values_from = .data$z) %>%
      dplyr::select(-.data$locus, -.data$allele) %>%
      as.matrix() %>% t()
  })
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # get distance between prediction sites
  grid_coords <- project$maps$grid$centroids
  dist_22 <- get_GC_distance(grid_coords$longitude, grid_coords$latitude) %>%
    as.matrix()
  
  # get distance between sampling sites and prediction sites
  site_coords <- project$data$raw$site_data
  dist_12 <- apply(grid_coords, 1, function(y) {
    lonlat_to_bearing(site_coords$longitude, site_coords$latitude, y[1], y[2])$gc_dist
  })
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, reps, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create misc list
  args_misc <- list(pb_markdown = pb_markdown)
  
  # initialise list for storing results
  project$maps$predictions <- replicate(length(loci), NULL)
  names(project$maps$predictions) <- sprintf("locus_%s", loci)
  
  # loop through loci
  for (i in seq_along(loci)) {
    message(sprintf("\nlocus %s of %s", i, length(loci)))
    
    # draw from predictive distribution via efficient C++ function
    output_raw <- predict_map_cpp(data_list[[i]], mcmc_sample, dist_11, dist_12,
                                  dist_22, project$model$parameters, inner_reps,
                                  args_progress, args_functions, args_misc)
    
    #return(output_raw)
    
    # get raw output into single array
    sim_array <- abind::abind(output_raw$ret, along = 2)
    
    # get mean and standard deviation over sims
    project$maps$predictions[[i]]$mean <-  apply(sim_array, c(1, 3), mean)
    project$maps$predictions[[i]]$sd <-  apply(sim_array, c(1, 3), sd)
    
    # get quantiles
    if (!is.null(quantiles)) {
      if (length(quantiles) == 1) {
        sim_quants <- apply(sim_array, c(1, 3), quantile, probs = quantiles) %>%
          list()
      } else {
        sim_quants <- apply(sim_array, c(1, 3), quantile, probs = quantiles) %>%
          purrr::array_tree(margin = 1)
      }
      names(sim_quants) <- sprintf("%s%%", round(quantiles * 100, digits = 1))
      project$maps$predictions[[i]]$quantiles <-  sim_quants
    }
    
    # get exceedance
    if (!is.null(exceedance)) {
      sim_exceedance <- mapply(function(x) {
        apply(sim_array > x, c(1, 3), mean)
      }, exceedance, SIMPLIFY = FALSE)
      names(sim_exceedance) <- sprintf("%s%%", round(exceedance * 100, digits = 1))
      project$maps$predictions[[i]]$exceedance <-  sim_exceedance
    }
    
  }  # end loop over loci
  
  return(project)
}

#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param reps TODO
#' @param inner_reps TODO
#' @param quantiles TODO
#' @param silent TODO
#' @param pb_markdown TODO
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats quantile sd density
#' @export

predict_pairwise <- function(project, reps = 2, inner_reps = 10,
                             quantiles = c(0.025, 0.5, 0.975),
                             silent = FALSE, pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_single_pos_int(inner_reps, zero_allowed = FALSE)
  if (!is.null(quantiles)) {
    assert_vector_bounded(quantiles)
  }
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # get observed Gst values
  obs_Gst <- as.matrix(project$data$pairwise_measures$Gst)
  n_site <- nrow(obs_Gst)
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$lambda, .data$nu) %>%
    dplyr::sample_n(reps) %>%
    as.list()
  
  # transform frequencies to continuous scale
  data_df <- project$data$raw$genetic_data %>%
    transform_p_to_z()
  
  # get z values split by locus and grouped into matrix by allele
  data_list <- lapply(split(data_df, data_df$locus), function(x) {
    tidyr::pivot_wider(x, names_from = .data$site_ID, values_from = .data$z) %>%
      dplyr::select(-.data$locus, -.data$allele) %>%
      as.matrix() %>% t()
  })
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, reps, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create misc list
  args_misc <- list(silent = silent,
                    pb_markdown = pb_markdown)
  
  # draw from predictive distribution via efficient C++ function
  output_raw <- predict_pairwise_Gst_cpp(data_list, mcmc_sample, dist_11,
                                         project$model$parameters, inner_reps,
                                         args_progress, args_functions, args_misc)
  
  #return(output_raw)
  
  # get raw output into single array
  sim_array <- abind::abind(output_raw$ret, along = 3)
  
  # get mean and sd of distribution
  if (!silent) {
    message("Calculating mean and SD of null distribution")
  }
  sim_mean <- sim_sd <- matrix(NA, n_site, n_site)
  for (i in 1:(n_site - 1)) {
    for (j in (i + 1):n_site) {
      sim_mean[i,j] <- mean(sim_array[i,j,])
      sim_sd[i,j] <- sd(sim_array[i,j,])
    }
  }
  
  # get quantiles of distribution
  sim_quantiles <- NULL
  if (!is.null(quantiles)) {
    if (!silent) {
      message("Calculating quantiles of null distribution")
    }
    sim_quantiles <- array(NA, dim = c(n_site, n_site, length(quantiles)))
    for (i in 1:(n_site - 1)) {
      for (j in (i + 1):n_site) {
        sim_quantiles[i,j,] <- quantile(sim_array[i,j,], probs = quantiles)
      }
    }
  }
  
  # make into list
  sim_quantiles <- mapply(function(i) sim_quantiles[,,i], seq_along(quantiles), SIMPLIFY = FALSE)
  names(sim_quantiles) <- sprintf("Q%s", round(quantiles * 100, digits = 1))
  
  # get ranking
  if (!silent) {
    message("Calculating ranking of observed values")
  }
  sim_rank <- matrix(NA, n_site, n_site)
  for (i in 1:(n_site - 1)) {
    for (j in (i + 1):n_site) {
      sim_rank[i,j] <- mean(obs_Gst[i,j] > sim_array[i,j,])
    }
  }
  
  # get empirical p-value via kernel density estimation
  if (!silent) {
    message("Calculating empirical p-values of observed values")
  }
  sim_p <- matrix(NA, n_site, n_site)
  for (i in 1:(n_site - 1)) {
    for (j in (i + 1):n_site) {
      bw <- density(sim_array[i,j,])$bw
      p_raw <- mean(pnorm(obs_Gst[i,j], mean = sim_array[i,j,], sd = bw))
      sim_p[i,j] <- 1 - 2 * abs(p_raw - 0.5)
    }
  }
  
  # load into project
  project$pairwise_predictions$Gst <- list(mean = sim_mean,
                                           sd = sim_sd,
                                           quantiles = sim_quantiles,
                                           ranking = sim_rank,
                                           empirical_p = sim_p)
  
  return(project)
}

#------------------------------------------------
#' @title Assign edges for GeoMAPI analysis
#'
#' @description Given an GeneScapeR project with a grid already loaded (see
#'   \code{?create_hex_grid}), determine which pairwise edges in the data
#'   intersect each grid cell. Assumes an elliptical projection along each edge
#'   with the start- and end-points becoming the two foci of the ellipse, and
#'   with the eccentricity defined by the user. See Details for suggestions on
#'   how to choose this parameter.
#'
#' @details TODO
#' 
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param eccentricity eccentricity of ellipses, defined as half the distance
#'   between foci divided by the semi-major axis. \eqn{e = sqrt{1 - b^2/a^2}},
#'   where \eqn{e} is the eccentricity, \eqn{a} is the length of the semi-major
#'   axis, and \eqn{b} is the length of the semi-minor axis. Eccentricity ranges
#'   between 0 (perfect circle) and 1 (straight line between foci). An
#'   eccentricity of 0 is not allowed in this case because it would result in an
#'   infinitely large circle.
#' @param max_dist edges shorter than this length are discarded prior to
#'   assignment. Can be used to focus on short distance signals.
#' @param silent if \code{TRUE} then no output is produced during function
#'   evaluation.
#' @param pb_markdown whether to run progress bars in markdown mode, in which
#'   case they are updated once at the end to avoid large amounts of output.
#'
#' @importFrom stats as.dist
#' @importFrom utils txtProgressBar
#' @export

GeoMAPI_assign_edges <- function(project,
                                 eccentricity = 0.9,
                                 max_dist = Inf,
                                 silent = FALSE,
                                 pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_bounded(eccentricity, inclusive_left = FALSE, inclusive_right = TRUE)
  assert_single_numeric(max_dist)
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # deal with infinite max_dist
  if (!is.finite(max_dist)) {
    max_dist <- max(dist_11) +1e2
  }
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  args_progress <- list()
  if (!silent) {
    n_poly <- length(project$maps$grid$polygons)
    pb <- txtProgressBar(0, n_poly, initial = NA, style = 3)
    args_progress <- list(pb = pb)
  }
  
  # create argument list
  args <- list(node_lon = project$data$raw$site_data$longitude,
               node_lat = project$data$raw$site_data$latitude,
               centroid_lon = project$maps$grid$centroids$longitude,
               centroid_lat = project$maps$grid$centroids$latitude,
               width = project$maps$grid$parameters$hex_width,
               eccentricity = eccentricity,
               max_dist = max_dist,
               silent = silent,
               pb_markdown = pb_markdown)
  
  # assign edges via efficient C++ function
  output_raw <- GeoMAPI_assign_edges_cpp(args, args_functions, args_progress, dist_11)
  
  # calculate coverage
  coverage <- mapply(length, output_raw$edge_assignment)
  
  # save results to project
  project$GeoMAPI$parameters <- list(eccentricity = eccentricity)
  project$GeoMAPI$edge_assignment <- output_raw$edge_assignment
  project$GeoMAPI$coverage <- coverage
  
  return(project)
}

#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param reps TODO
#' @param inner_reps TODO
#' @param silent TODO
#' @param pb_markdown TODO
#'
#' @importFrom utils txtProgressBar
#' @export

GeoMAPI_analysis <- function(project, reps = 2, inner_reps = 10,
                             silent = FALSE, pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_single_pos_int(inner_reps, zero_allowed = FALSE)
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # get basic dimensions
  n_site <- nrow(project$data$raw$site_data)
  n_cells <- project$GeoMAPI$edge_assignment %>% length()
  total_reps <- reps * inner_reps
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$lambda, .data$nu) %>%
    dplyr::sample_n(reps) %>%
    as.list()
  
  # transform frequencies to continuous scale
  data_df <- project$data$raw$genetic_data %>%
    transform_p_to_z()
  
  # get z values split by locus and grouped into matrix by allele
  data_list <- lapply(split(data_df, data_df$locus), function(x) {
    tidyr::pivot_wider(x, names_from = .data$site_ID, values_from = .data$z) %>%
      dplyr::select(-.data$locus, -.data$allele) %>%
      as.matrix() %>% t()
  })
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, reps, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create misc list
  args_misc <- list(silent = silent,
                    pb_markdown = pb_markdown)
  
  # draw from predictive distribution via efficient C++ function
  output_raw <- predict_pairwise_Gst_cpp(data_list, mcmc_sample, dist_11,
                                         project$model$parameters, inner_reps,
                                         args_progress, args_functions, args_misc)
  
  # get raw output into single array
  sim_array <- abind::abind(output_raw$ret, along = 3)
  
  # create matrix relating edges to nodes
  mat_node_edge <- cbind(node_1 = rep(1:(n_site - 1), times = (n_site - 1):1),
                         node_2 = mapply(function(x) x:n_site, 2:n_site) %>% unlist())
  
  # get observed Gst matrix
  obs_Gst <- as.matrix(project$data$pairwise_measures$Gst)
  
  if (!silent) {
    message("Calculating z_scores")
  }
  
  # calculate z_score for all grid cells. Nb, this was found to be no faster in
  # C++
  z_score <- rep(NA, n_cells)
  for (i in 1:n_cells) {
    w_e <- project$GeoMAPI$edge_assignment[[i]]
    n_e <- length(w_e)
    if (n_e == 0) {
      next()
    }
    w_n <- mat_node_edge[w_e,, drop = FALSE]
    obs_stat <- mean(obs_Gst[w_n])
    sim_stat <- 0
    for (j in seq_len(n_e)) {
      sim_stat <- sim_stat + sim_array[w_n[j,1], w_n[j,2],]
    }
    sim_stat <- sim_stat / n_e
    z_score[i] <- (obs_stat - mean(sim_stat)) / sd(sim_stat)
  }
  
  # save to project
  project$GeoMAPI$z_score <- z_score
  
  return(project)
}

#------------------------------------------------
#' @title Get GeoMAPI significant cells
#'
#' @description Given a completed GeoMAPI analysis, identify which cells are
#'   significant outliers.
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param test_tail whether to calculate empirical p-values using a
#'   one-sided test (\code{test_tail = "left"} or \code{test_tail =
#'   "right"}) or a two-sided test (\code{test_tail = "both"}).
#' @param FDR the false discovery rate, i.e. the probability that a cell
#'   identified as significant is actually a false positive.
#' @param min_coverage minimum coverage (number of edges assigned to a cell)
#'   for it to be included in the final result.
#'
#' @importFrom stats pnorm
#' @export

GeoMAPI_get_significant <- function(project,
                                    test_tail = "both",
                                    FDR = 0.05,
                                    min_coverage = 10) {
  
  # check project
  assert_class(project, "genescaper_project")
  assert_single_string(test_tail)
  assert_in(test_tail, c("left", "right", "both"))
  assert_single_bounded(FDR)
  assert_single_pos_int(min_coverage, zero_allowed = TRUE)
  
  # get results into dataframe
  df_res <- data.frame(cell = seq_along(project$maps$grid$polygons),
                       z_score = project$GeoMAPI$z_score,
                       coverage = project$GeoMAPI$coverage)
  
  # subset based on coverage
  if (!any(df_res$coverage >= min_coverage)) {
    stop("no cells wth sufficient coverage when calculating significance")
  }
  df_res <- dplyr::filter(df_res, .data$coverage >= min_coverage)
  
  # calculate p-values
  if (test_tail == "left") {
    df_res$p <- pnorm(df_res$z_score) 
  } else if (test_tail == "right") {
    df_res$p <- pnorm(df_res$z_score, lower.tail = FALSE) 
  } else if (test_tail == "both") {
    df_res$p <- 2*pnorm(-abs(df_res$z_score))
  }
  
  # sort in order of increasing p
  df_res <- df_res[order(df_res$p),]
  
  # Bejamini and Yekutieli (2001) method for identifying significant results
  # while fixing the false descovery rate
  df_res$BY <- FDR * seq_along(df_res$p) / nrow(df_res)
  which_lower <- which_upper <- integer()
  if (any(df_res$p <= df_res$BY, na.rm = TRUE)) {
    
    w <- which(df_res$p <= df_res$BY)
    which_upper <- df_res$cell[w][df_res$z_score[w] > 0]
    which_lower <- df_res$cell[w][df_res$z_score[w] <= 0]
  }
  
  # add to project
  project$GeoMAPI$significance <- list(upper = which_upper,
                                       lower = which_lower)
  
  return(project)
}
