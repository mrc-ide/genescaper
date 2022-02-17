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
  assert_in(c("site_ID", "latitude", "longitude"), names(site_data),
            message = "site_data must have columns {site_ID, latitude, longitude}")
  assert_vector_numeric(site_data$latitude)
  assert_bounded(site_data$latitude, left = -90, right = 90)
  assert_vector_numeric(site_data$longitude)
  assert_bounded(site_data$longitude, left = -180, right = 180)
  
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
  site_dist <- get_GC_distance(project$data$raw$site_data$longitude,
                               project$data$raw$site_data$latitude) %>%
    as.matrix()
  
  # get adjusted allele frequencies by applying stick breaking correction. For
  # each allele, p_stick is calculated by dividing the allele frequency by the
  # amount of unit interval that remains once previous frequencies have been
  # taken into account. This creates (J - 1) ostensibly independent frequencies,
  # where J is the number of alleles. The final p_stick always equals 1 and so
  # this allele can be dropped (in other words, there are J - 1 degrees of
  # freedom and so we are reducing from J dependent frequencies to J-1
  # independent frequencies).
  # z values are calculated as logit(p_stick). Finally, z values are split into
  # a list over all locus&allele combos. Each combo has an independent mean and
  # variance under the model and so we can treat these as equivalent replicates.
  z_list <- project$data$raw$genetic_data %>%
    dplyr::group_by(.data$site_ID, .data$locus) %>%
    dplyr::summarise(J = length(.data$allele),
                     allele = .data$allele[-.data$J],
                     stick_remaining = 1 - cumsum(.data$freq[-.data$J]) + .data$freq[-.data$J],
                     p_stick = .data$freq[-.data$J] / .data$stick_remaining,
                     z = log(.data$p_stick) - log(1 - .data$p_stick)) %>%
    dplyr::ungroup() %>%
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
    assert_in(c("long", "lat"), names(border_coords))
    assert_vector_numeric(border_coords$longitude)
    assert_vector_numeric(border_coords$latitude)
    assert_bounded(border_coords$longitude, left = -180, right = 180)
    assert_bounded(border_coords$latitude, left = -90, right = 90)
  }
  
  message("Creating hex map")
  
  # calculate default hex size from data
  if (is.null(hex_width)) {
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
  names(hex_pts_df) <- c("long", "lat")
  
  message(sprintf("%s hexagons created", nhex))
  
  # add to project
  project$maps <- list(grid = list(parameters = list(hex_width = hex_width),
                                   centroids = hex_pts_df,
                                   polygons = hex_polys))
  
  return(project)
}
