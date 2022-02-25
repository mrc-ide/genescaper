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
  
  # get allele frequencies into list over loci then alleles
  freq_array <- project$data$raw$genetic_data %>%
    group_by(locus) %>%
    group_split() %>%
    lapply(function (x) {
      x %>% group_by(allele) %>%
        group_split() %>%
        lapply(function (x) x$freq)
    })
  
  # pass to efficient C++ function
  output_raw <- wrangle_pairwise_Gst_cpp(freq_array)
  
  # convert to distance matrix
  n_site <- nrow(project$data$raw$site_data)
  Gst_mat <- matrix(0, n_site, n_site)
  Gst_mat[lower.tri(Gst_mat)] <- output_raw$Gst
  project$data$pairwise_measures$Gst <- as.dist(Gst_mat)
  
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
  site_dist <- as.matrix(myproj$data$pairwise_measures$distance)
  
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
  
  # get z values split by locus and allele
  data_list <- lapply(split(data_df, data_df$locus, drop = TRUE),
                      function(x) split(x[["z"]], x[['allele']], drop = TRUE))
  
  # get distance between sampling sites
  site_coords <- project$data$raw$site_data %>%
    dplyr::select(.data$longitude, .data$latitude)
  dist_11 <- get_GC_distance(site_coords$longitude, site_coords$latitude) %>%
    as.matrix()
  
  # get distance between prediction sites
  grid_coords <- project$maps$grid$centroids
  dist_22 <- get_GC_distance(grid_coords$longitude, grid_coords$latitude) %>%
    as.matrix()
  
  # get distance between sampling sites and prediction sites
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
    
    # get raw output into array
    sim_array <- mapply(function(x) {
      matrix(unlist(x), ncol = length(x))
    }, output_raw$ret, SIMPLIFY = "array")
    
    # get mean and standard deviation over sims
    project$maps$predictions[[i]]$mean <-  apply(sim_array, c(1, 2), mean)
    project$maps$predictions[[i]]$sd <-  apply(sim_array, c(1, 2), sd)
    
    # get quantiles
    if (!is.null(quantiles)) {
      if (length(quantiles) == 1) {
        sim_quants <- apply(sim_array, c(1, 2), quantile, probs = quantiles) %>%
          list()
      } else {
        sim_quants <- apply(sim_array, c(1, 2), quantile, probs = quantiles) %>%
          purrr::array_tree(margin = 1)
      }
      names(sim_quants) <- sprintf("%s%%", round(quantiles * 100, digits = 1))
      project$maps$predictions[[i]]$quantiles <-  sim_quants
    }
    
    # get exceedance
    if (!is.null(exceedance)) {
      sim_exceedance <- mapply(function(x) {
        apply(sim_array > x, c(1, 2), mean)
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
#' @param pb_markdown TODO
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats quantile sd
#' @export

predict_pairwise <- function(project, reps = 2, inner_reps = 10,
                        quantiles = c(0.025, 0.5, 0.975), pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_single_pos_int(inner_reps, zero_allowed = FALSE)
  if (!is.null(quantiles)) {
    assert_vector_bounded(quantiles)
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
    transform_p_to_z()
  
  # get z values split by locus and allele
  data_list <- lapply(split(data_df, data_df$locus, drop = TRUE),
                      function(x) split(x[["z"]], x[['allele']], drop = TRUE))
  
  # get distance between sampling sites
  site_coords <- project$data$raw$site_data %>%
    dplyr::select(.data$longitude, .data$latitude)
  dist_11 <- get_GC_distance(site_coords$longitude, site_coords$latitude) %>%
    as.matrix()
  
  # get distance between prediction sites
  grid_coords <- project$maps$grid$centroids
  dist_22 <- get_GC_distance(grid_coords$longitude, grid_coords$latitude) %>%
    as.matrix()
  
  # get distance between sampling sites and prediction sites
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
    
    # get raw output into array
    sim_array <- mapply(function(x) {
      matrix(unlist(x), ncol = length(x))
    }, output_raw$ret, SIMPLIFY = "array")
    
    # get mean and standard deviation over sims
    project$maps$predictions[[i]]$mean <-  apply(sim_array, c(1, 2), mean)
    project$maps$predictions[[i]]$sd <-  apply(sim_array, c(1, 2), sd)
    
    # get quantiles
    if (!is.null(quantiles)) {
      if (length(quantiles) == 1) {
        sim_quants <- apply(sim_array, c(1, 2), quantile, probs = quantiles) %>%
          list()
      } else {
        sim_quants <- apply(sim_array, c(1, 2), quantile, probs = quantiles) %>%
          purrr::array_tree(margin = 1)
      }
      names(sim_quants) <- sprintf("%s%%", round(quantiles * 100, digits = 1))
      project$maps$predictions[[i]]$quantiles <-  sim_quants
    }
    
    # get exceedance
    if (!is.null(exceedance)) {
      sim_exceedance <- mapply(function(x) {
        apply(sim_array > x, c(1, 2), mean)
      }, exceedance, SIMPLIFY = FALSE)
      names(sim_exceedance) <- sprintf("%s%%", round(exceedance * 100, digits = 1))
      project$maps$predictions[[i]]$exceedance <-  sim_exceedance
    }
    
  }  # end loop over loci
  
  return(project)
}
