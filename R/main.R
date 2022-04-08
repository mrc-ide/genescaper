#------------------------------------------------
#' @title Load data into project
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param site_data TODO
#' @param genetic_data TODO
#' @param epsilon TODO
#'
#' @importFrom rlang .data
#' @export

bind_data <- function(project, site_data, genetic_data, epsilon = 0.5) {
  
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
  mssg1 <- paste0("genetic_data must have columns {site_ID, locus, allele}, along with either {freq}",
                  " for frequency data, or {count} for count data.")
  assert_in(c("site_ID", "locus", "allele"), names(genetic_data), message = mssg1)
  if (!("freq" %in% names(genetic_data)) & !("count" %in% names(genetic_data))) {
    stop(mssg1)
  }
  if (("freq" %in% names(genetic_data)) & ("count" %in% names(genetic_data))) {
    stop("genetic_data must contain columns {freq} OR {count}, but not both.")
  }
  assert_vector_pos_int(genetic_data$locus)
  assert_vector_pos_int(genetic_data$allele)
  genetic_type = ifelse("freq" %in% names(genetic_data), "freq", "count")
  if (genetic_type == "freq") {
    assert_vector_pos(genetic_data$freq)
  }
  if (genetic_type == "count") {
    assert_vector_pos_int(genetic_data$count, zero_allowed = TRUE)
  }
  
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
  
  # check that frequencies sum to 1
  if (genetic_type == "freq") {
    freq_sum <- genetic_data %>%
      dplyr::group_by(.data$site_ID, .data$locus) %>%
      dplyr::summarise(freq_sum = sum(.data$freq))
    if (!isTRUE(all.equal(freq_sum$freq_sum, rep(1, nrow(freq_sum))))) {
      stop("allele frequencies must sum to 1 at every site and every locus")
    }
  }
  
  # check that counts sum to non-zero value
  if (genetic_type == "count") {
    count_sum <- genetic_data %>%
      dplyr::group_by(.data$site_ID, .data$locus) %>%
      dplyr::summarise(count_sum = sum(.data$count))
    if (any(count_sum$count_sum == 0)) {
      stop("must be at least one non-zero count over all alleles at a locus")
    }
  }
  
  # process data to deal with allele frequencies of exactly 0 or 1
  if (genetic_type == "freq") {
    
    # check that frequencies are not fixed everywhere at any locus
    any_unfixed <- genetic_data %>%
      dplyr::group_by(.data$locus, .data$allele) %>%
      dplyr::summarise(any_unfixed = any((.data$freq != 0) & (.data$freq != 1)))
    if (!all(any_unfixed$any_unfixed)) {
      stop("frequencies cannot be fixed (frequency of exactly 0 or 1) at every spatial location")
    }
    
    # get epsilon value for each locus-allele combination
    df_epsilon <- genetic_data %>%
      dplyr::group_by(.data$locus, .data$allele) %>%
      dplyr::summarise(min_freq = min(.data$freq[(.data$freq > 0) & (.data$freq < 1)]),
                       max_freq = max(.data$freq[(.data$freq > 0) & (.data$freq < 1)]),
                       epsilon = min(c(.data$min_freq, 1 - .data$max_freq))) %>%
      dplyr::select("locus", "allele", "epsilon")
    
    # apply epsilon correction
    genetic_data <- genetic_data %>%
      dplyr::left_join(df_epsilon, by = c("locus", "allele")) %>%
      dplyr::mutate(freq_raw = .data$freq,
                    freq = ifelse(.data$freq < epsilon, epsilon, ifelse(.data$freq > 1 - epsilon, 1 - epsilon, .data$freq))) %>%
      dplyr::ungroup() %>%
      dplyr::select("site_ID", "locus", "allele", "freq_raw", "freq")
    
  } else if (genetic_type == "count") {
    genetic_data <- genetic_data %>%
      dplyr::mutate(count_adjusted = .data$count + epsilon) %>%
      dplyr::group_by(.data$site_ID, .data$locus) %>%
      dplyr::summarise(allele = .data$allele,
                       count = .data$count,
                       freq = .data$count_adjusted / sum(.data$count_adjusted))
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
  
  # get allele frequencies into list over loci then matrix over demes and
  # alleles
  freq_list <- project$data$raw$genetic_data %>%
    dplyr::select(.data$site_ID, .data$locus, .data$allele, .data$freq) %>%
    dplyr::group_by(.data$locus) %>%
    dplyr::group_split() %>%
    lapply(function (x) {
      x %>% tidyr::pivot_wider(names_from = .data$allele, values_from = .data$freq) %>%
        dplyr::select(-.data$site_ID, -.data$locus) %>%
        as.matrix()
    })
  
  # calculate mean Gst over loci
  Gst <- 0
  for (i in seq_along(freq_list)) {
    Gst <- Gst + calc_pairwise_Gst(freq_list[[i]])
  }
  Gst <- Gst / length(freq_list)
  
  # get into matrix
  Gst_mat <- pairwise_to_mat(Gst)
  
  # load into project
  project$data$pairwise_measures$Gst <- Gst_mat
  
  return(project)
}

#------------------------------------------------
#' @title Calculate pairwise Jost's D between sites
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#'
#' @export

get_pairwise_D <- function(project) {
  
  # get allele frequencies into list over loci then matrix over demes and
  # alleles
  freq_list <- project$data$raw$genetic_data %>%
    dplyr::select(.data$site_ID, .data$locus, .data$allele, .data$freq) %>%
    dplyr::group_by(.data$locus) %>%
    dplyr::group_split() %>%
    lapply(function (x) {
      x %>% tidyr::pivot_wider(names_from = .data$allele, values_from = .data$freq) %>%
        dplyr::select(-.data$site_ID, -.data$locus) %>%
        as.matrix()
    })
  
  # calculate mean Jost's D over loci
  D <- 0
  for (i in seq_along(freq_list)) {
    D <- D + calc_pairwise_D(freq_list[[i]])
  }
  D <- D / length(freq_list)
  
  # get into matrix
  D_mat <- pairwise_to_mat(D)
  
  # load into project
  project$data$pairwise_measures$D <- D_mat
  
  return(project)
}

#------------------------------------------------
#' @title Define parameters of the spatial correlation model
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param nu_shape1 TODO
#' @param nu_shape2 TODO
#' @param lambda_shape TODO
#' @param lambda_rate TODO
#'
#' @export

define_model <- function(project, nu_shape1 = 1.0, nu_shape2 = 1.0, lambda_shape = 1.0, lambda_rate = NULL) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos(nu_shape1, zero_allowed = FALSE)
  assert_single_pos(nu_shape2, zero_allowed = FALSE)
  assert_single_pos(lambda_shape, zero_allowed = FALSE)
  if (!is.null(lambda_rate)) {
    assert_single_pos(lambda_rate, zero_allowed = FALSE)
  }
  
  # set default lambda_rate such that expected lambda is half maximum distance
  # in data
  if (is.null(lambda_rate)) {
    max_dist <- max(project$data$pairwise_measures$distance)
    lambda_rate <- lambda_shape / (max_dist / 2)
  }
  
  # store within project
  project$model$parameters <- list(nu_shape1 = nu_shape1,
                                   nu_shape2 = nu_shape2,
                                   lambda_shape = lambda_shape,
                                   lambda_rate = lambda_rate)
  
  return(project)
}

#------------------------------------------------
#' @title Run MCMC
#'
#' @description TODO
#'
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param chains TODO
#' @param ... additional parameters that will be passed to
#'   \code{drjacoby::run_mcmc()}.
#'
#' @export

run_mcmc <- function(project, true_sigsq, chains, ...) {
  
  # avoid "no visible binding" warning
  loglike <- NULL
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(chains, zero_allowed = FALSE)
  
  # check that both data and model have been defined
  assert_non_null(project$dat$raw)
  assert_non_null(project$model$parameters)
  
  # get distance between all sites
  site_dist <- project$data$pairwise_measures$distance
  
  # transform allele frequencies to z values
  z_df <- project$data$raw$genetic_data %>%
    transform_p_to_z()
  
  # split z-values into a list over all locus-allele combos
  z_list <-  z_df %>%
    dplyr::select(.data$locus, .data$allele, .data$z) %>%
    dplyr::group_by(.data$locus, .data$allele) %>%
    dplyr::group_split() %>%
    lapply(function (x) x$z)
  names(z_list) <- seq_along(z_list)
  
  # define parameters dataframe
  param_type <- 3
  if (param_type == 1) {
    df_params <- drjacoby::define_params(name = "nu", min = 0, max = 1, init = rep(1e-4, chains),
                                         name = "log_lambda", min = -Inf, max = Inf, init = rep(log(mean(site_dist)), chains),
                                         name = "omega", min = 1, max = 3.0, init = runif(chains, 1.1, 1.5),
                                         name = "gamma", min = 0.1, max = 10, init = runif(chains, 0.2, 9))
  } else if (param_type == 2) {
    nu_fix <- nu_true
    lambda_fix <- lambda_true
    omega_fix <- omega_true
    gamma_fix <- gamma_true
    df_params <- drjacoby::define_params(name = "nu", min = 0, max = 1, init = rep(nu_fix, chains),
                                         name = "log_lambda", min = log(lambda_fix), max = log(lambda_fix), init = rep(log(lambda_fix), chains),
                                         name = "omega", min = omega_fix, max = omega_fix, init = rep(omega_fix, chains),
                                         name = "gamma", min = gamma_fix, max = gamma_fix, init = rep(gamma_fix, chains))
  } else if (param_type == 3) {
    nu_fix <- nu_true
    lambda_fix <- lambda_true
    omega_fix <- omega_true
    gamma_fix <- 1.0
    df_params <- drjacoby::define_params(name = "nu", min = 0, max = 1, init = rep(nu_fix, chains),
                                         name = "log_lambda", min = -Inf, max = Inf, init = rep(log(lambda_fix), chains),
                                         name = "omega", min = 1, max = 3, init = rep(omega_fix, chains),
                                         name = "gamma", min = gamma_fix, max = gamma_fix, init = rep(gamma_fix, chains))
  }
  
  # define misc list
  misc_list <- append(project$model$parameters,
                      list(site_dist = site_dist,
                           n_site = nrow(site_dist)))
  
  # TODO - remove
  misc_list$true_sigsq <- true_sigsq
  
  
  # source C++ likelihood and prior functions
  Rcpp::sourceCpp(system.file("extdata/GRF_model.cpp", package = 'genescaper', mustWork = TRUE))
  
  # check that all initial values create valid likelihoods
  ll_init <- rep(NA, chains)
  for (i in 1:chains) {
    param_vec <- c("nu" = df_params$init[[1]][i],
                   "log_lambda" = df_params$init[[2]][i],
                   "omega" = df_params$init[[3]][i],
                   "gamma" = df_params$init[[4]][i])
    
    ll_init[i] <- loglike(params = param_vec, data = z_list, misc = misc_list)
  }
  if (any(ll_init < -1e+300)) {
    stop("initial parameter values produce log-likelihoods outside reasonable range")
  }
  
  # run MCMC
  mcmc <- drjacoby::run_mcmc(data = z_list,
                             df_params = df_params,
                             misc = misc_list,
                             loglike = "loglike",
                             logprior = "logprior",
                             chains = chains,
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
  
  # get adjacency list
  hex_adj <- st_touches(hex_polys)
  
  message(sprintf("%s hexagons created", nhex))
  
  # add to project
  project$maps <- list(grid = list(parameters = list(hex_width = hex_width),
                                   centroids = hex_pts_df,
                                   polygons = hex_polys,
                                   adjacency = hex_adj))
  
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
#' @param silent TODO
#' @param pb_markdown TODO
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats quantile sd
#' @export

predict_map <- function(project, loci, reps = 2, inner_reps = 10,
                        quantiles = c(0.025, 0.5, 0.975), exceedance = NULL,
                        silent = FALSE, pb_markdown = FALSE) {
  
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
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$nu, .data$lambda, .data$omega, .data$gamma) %>%
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
  args_misc <- list(silent = silent,
                    pb_markdown = pb_markdown)
  
  # initialise list for storing results
  project$maps$predictions <- replicate(length(loci), NULL)
  names(project$maps$predictions) <- sprintf("locus_%s", loci)
  
  # loop through loci
  for (i in seq_along(loci)) {
    if (!silent) {
      message(sprintf("\nlocus %s of %s", i, length(loci)))
    }
    
    # draw from predictive distribution via efficient C++ function
    sim_array <- predict_map_cpp(data_list[[i]], mcmc_sample, dist_11, dist_12,
                                  dist_22, project$model$parameters, inner_reps,
                                  args_progress, args_functions, args_misc) %>%
      abind::abind(along = 2)
    
    #return(sim_array)
    
    if (!silent) {
      if (i == 1) {
        message("processing")
      } else {
        message("\nprocessing")
      }
    }
    
    # get mean and standard deviation over sims
    project$maps$predictions[[i]]$mean <-  apply(sim_array, c(1, 3), mean)
    project$maps$predictions[[i]]$sd <-  apply(sim_array, c(1, 3), sd)
    
    # get quantiles
    if (!is.null(quantiles)) {
      sim_quants <- apply(sim_array, c(1, 3), quantile, probs = quantiles)
      if (length(quantiles) == 1) {
        sim_quants <- list(sim_quants)
      } else {
        sim_quants <- purrr::array_tree(sim_quants, margin = 1)
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
  
  # simulate Gst from null model
  Gst_sim <- stat_sim(project, reps, inner_reps, silent, pb_markdown)
  
  if (!silent) {
    message("Calculating summaries")
  }
  
  # get observed Gst values
  Gst_obs <- project$data$pairwise_measures$Gst
  n_site <- nrow(Gst_obs)
  
  # get mean and standard deviation over sims
  sim_mean <- pairwise_to_mat(rowMeans(Gst_sim))
  sim_sd <- pairwise_to_mat(apply(Gst_sim, 1, sd))
  
  # get quantiles
  if (!is.null(quantiles)) {
    sim_quants <- list()
    for (i in seq_along(quantiles)) {
      sim_quants[[i]] <- pairwise_to_mat(apply(Gst_sim, 1, quantile, probs = quantiles[i]))
    }
    names(sim_quants) <- sprintf("%s%%", round(quantiles * 100, digits = 1))
  }
  
  # get percentile rank of observed data
  Gst_obs_vec <- Gst_obs[lower.tri(Gst_obs)]
  Gst_obs_mat <- matrix(Gst_obs_vec, nrow = length(Gst_obs_vec), ncol = reps * inner_reps)
  sim_percentile_rank <- (rowSums(Gst_obs_mat > Gst_sim) + 1) / (reps * inner_reps + 1) * 100
  
  # load into project
  project$pairwise_predictions$Gst <- list(mean = sim_mean,
                                           sd = sim_sd,
                                           quantiles = sim_quants,
                                           percentile_rank = sim_percentile_rank)
  
  return(project)
}

# -----------------------------------
# simulates pairwise stats by drawing from null model
#' @noRd
stat_sim <- function(project, reps, inner_reps, silent, pb_markdown) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_single_pos_int(inner_reps, zero_allowed = FALSE)
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$nu, .data$lambda, .data$omega, .data$gamma) %>%
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
  loci <- seq_along(data_list)
  n_loci <- length(loci)
  alleles <- mapply(ncol, data_list)
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, n_loci, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create misc list
  args_misc <- list(silent = TRUE,
                    pb_markdown = TRUE)
  
  if (!silent) {
    message("Drawing from null distribution for each locus")
    update_progress(args_progress, "pb", 0, n_loci)
  }
  
  # loop through loci
  ret <- 0
  for (i in seq_along(loci)) {
    
    # draw from null distribution via efficient C++ function
    sim_array <- null_site_cpp(data_list[[i]], mcmc_sample, dist_11,
                               project$model$parameters, inner_reps,
                               args_progress, args_functions, args_misc) %>%
      abind::abind(along = 2)
    
    # add to running estimate
    ret <- ret + calc_pairwise_Gst(sim_array)
    
    if (!silent) {
      update_progress(args_progress, "pb", i, n_loci)
    }
  } # end loop through loci
  
  # divide through by loci
  ret <- ret / n_loci
  
  return(ret)
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
    max_dist <- max(dist_11) + 1e2
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
#' @title Suggest eccentricity parameter based on target coverage
#'
#' @description TODO
#'
#' @details TODO
#' 
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function.
#' @param target_coverage,target_proportion eccentricity will be chosen such
#'   that a proportion \code{target_proportion} of cells achieve a coverage of
#'   at least \code{target_coverage}.
#' @param n_iterations the number of iterations to run of the binary search
#'   method.
#' @param max_dist edges shorter than this length are discarded prior to
#'   assignment. Can be used to focus on short distance signals.
#' @param silent if \code{TRUE} then no output is produced during function
#'   evaluation.
#' @param pb_markdown whether to run progress bars in markdown mode, in which
#'   case they are updated once at the end to avoid large amounts of output.
#'
#' @export

GeoMAPI_suggest_eccentricity <- function(project,
                                         target_coverage = 10,
                                         target_proportion = 0.9,
                                         n_iterations = 30,
                                         max_dist = Inf,
                                         silent = FALSE,
                                         pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(target_coverage, zero_allowed = FALSE)
  assert_single_bounded(target_proportion, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_pos_int(n_iterations, zero_allowed = FALSE)
  assert_single_numeric(max_dist)
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # get largest and smallest distance between sampling sites
  d_vec <- dist_11[upper.tri(dist_11)]
  d_min <- min(d_vec[d_vec > 0])
  d_max <- max(d_vec[d_vec > 0])
  
  # deal with infinite max_dist (cannot pass infinite values to C++)
  if (!is.finite(max_dist)) {
    max_dist <- d_max + 1e2
  }
  
  # get distance between cells, and get maximum distance
  cell_dist <- get_GC_distance(lon = project$maps$grid$centroids$longitude,
                               lat = project$maps$grid$centroids$latitude)
  c_max <- max(cell_dist)
  
  # use min/max distances to define a value of eccentricity that we know for
  # certain will result in all edges being assigned to all cells. The logic here
  # is to imagine a circle within an ellipse, with radius of the circle equal to
  # the semi-minor axis of the ellipse. Imagine that we choose the radius of the
  # circle to be the furthest distance between any two cells (c_max), such that
  # no matter where the circle is centred it will be assigned to all cells. This
  # tells us the desired semi-minor axis of the ellipse. Now imagine that the
  # ellipse is constructed between the *closest* two sampling locations. This
  # tells us the half-distance between foci. These two values together can be
  # used to derive the eccentricity.
  ecc_min <- 1 / sqrt((2 * c_max / d_min)^2 + 1)
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create argument list
  args <- list(node_lon = project$data$raw$site_data$longitude,
               node_lat = project$data$raw$site_data$latitude,
               centroid_lon = project$maps$grid$centroids$longitude,
               centroid_lat = project$maps$grid$centroids$latitude,
               width = project$maps$grid$parameters$hex_width,
               eccentricity = 1.0,
               max_dist = max_dist,
               silent = TRUE,
               pb_markdown = TRUE)
  
  # run coverage analysis once with eccentricity = 1.0, i.e. straight lines
  # between cells. Establish the proportion of cells with satisfactory coverage
  # in this situation. If it is greater than the desired target proportion then
  # there is no point performing a search, as eccentricity will simply keep
  # pushing up towards 1.
  output_raw <- GeoMAPI_assign_edges_cpp(args, args_functions, list(), dist_11)
  coverage <- mapply(length, output_raw$edge_assignment)
  actual_proportion <- mean(coverage > target_coverage)
  if (actual_proportion > target_proportion) {
    message(paste0("Even an eccentricity of 1.0 (i.e. perfect straight lines) achieves the target proportion. Think about:",
                   "\n - decreasing the max_dist (sharper resolution)",
                   "\n - increasing the target_proportion (greater number of reliable cells)",
                   "\n - increasing the target_coverage (more stringent requirement of cells)",
                   "\n - increasing the resolution of the grid (sharper resolution)",
                   "\n - setting the eccentricity manually"))
    return(1.0)
  }
  
  # set starting eccentricity bounds in search
  ecc_left <- ecc_min
  ecc_right <- 1.0
  ecc_prop <- ecc_min
  actual_proportion <- 1
  
  # create progress bar
  pb_search <- txtProgressBar(0, n_iterations, initial = NA, style = 3)
  
  # perform search
  message("Running binary search")
  for (i in 1:n_iterations) {
    
    # calculate next search value
    if (actual_proportion > target_proportion) {
      ecc_left <- ecc_prop
    } else {
      ecc_right <- ecc_prop
    }
    ecc_prop <- (ecc_left + ecc_right) / 2
    
    # update arguments with proposed eccentricity
    args$eccentricity <- ecc_prop
    
    # assign edges via efficient C++ function
    output_raw <- GeoMAPI_assign_edges_cpp(args, args_functions, list(), dist_11)
    
    # calculate coverage
    coverage <- mapply(length, output_raw$edge_assignment)
    
    # get proportion good coverage
    actual_proportion <- mean(coverage > target_coverage)
    
    # update progress bar
    if (!silent) {
      update_progress(list(pb = pb_search), "pb", i, n_iterations)
    }
  }
  
  # report result
  if (!silent) {
    message(sprintf("Eccentricity of %s means that %s%% of cells achieve a coverage of %s or higher",
                    ecc_prop, signif(actual_proportion * 100), target_coverage))
  }
  
  return(ecc_prop)
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
  n_cells <- project$GeoMAPI$edge_assignment %>% length()
  
  # get observed Gst
  Gst_obs_mat <- project$data$pairwise_measures$Gst
  Gst_obs <- Gst_obs_mat[lower.tri(Gst_obs_mat)]
  
  # simulate Gst from null model
  Gst_sim <- stat_sim(project, reps, inner_reps, silent, pb_markdown)
  
  if (!silent) {
    message("Calculating z_scores")
  }
  
  # normalise observed and simulated values
  Gst_mean <- rowMeans(Gst_sim)
  Gst_sd <- apply(Gst_sim, 1, sd)
  Gst_obs_norm <- (Gst_obs - Gst_mean) / Gst_sd
  Gst_sim_norm <- apply(Gst_sim, 2, function(x) {
    (x - Gst_mean) / Gst_sd
  })
  
  # calculate z_score for all grid cells
  z_score <- rep(NA, n_cells)
  for (i in 1:n_cells) {
    edges <- project$GeoMAPI$edge_assignment[[i]]
    n_edges <- length(edges)
    if (n_edges == 0) {
      next()
    }
    cell_obs <- mean(Gst_obs_norm[edges])
    cell_sim <- colMeans(Gst_sim_norm[edges,,drop = FALSE])
    z_score[i] <- (cell_obs - mean(cell_sim)) / sd(cell_sim)
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
#' @param silent TODO
#'
#' @importFrom stats pnorm
#' @export

GeoMAPI_get_significant <- function(project, test_tail = "both", FDR = 0.05,
                                    min_coverage = 10, silent = FALSE) {
  
  # check project
  assert_class(project, "genescaper_project")
  assert_single_string(test_tail)
  assert_in(test_tail, c("left", "right", "both"))
  assert_single_bounded(FDR)
  assert_single_pos_int(min_coverage, zero_allowed = TRUE)
  assert_single_logical(silent)
  
  # get results into dataframe
  df_res <- data.frame(cell = seq_along(project$maps$grid$polygons),
                       z_score = project$GeoMAPI$z_score,
                       coverage = project$GeoMAPI$coverage)
  
  # subset based on coverage
  if (!any(df_res$coverage >= min_coverage)) {
    stop("no cells wth sufficient coverage when calculating significance")
  }
  df_res <- dplyr::filter(df_res, .data$coverage >= min_coverage)
  
  if (!silent) {
    message("Calculating significance")
  }
  
  # calculate p-values
  if (test_tail == "left") {
    df_res$p <- pnorm(df_res$z_score) 
  } else if (test_tail == "right") {
    df_res$p <- pnorm(df_res$z_score, lower.tail = FALSE) 
  } else if (test_tail == "both") {
    df_res$p <- 2*pnorm(-abs(df_res$z_score))
  }
  
  # Bejamini and Yekutieli (2001) method for identifying significant results
  # while fixing the false descovery rate
  df_res$direction <- df_res$z_score
  which_signif <- Bejamini_Yekutieli(df_res, FDR)
  
  # add to project
  project$GeoMAPI$significance <- which_signif
  
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
#' @param measure TODO
#' @param patch_size TODO
#' @param quantiles TODO
#' @param silent TODO
#' @param pb_markdown TODO
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats quantile sd var
#' @export

Wombling <- function(project, loci = NULL, reps = 2, inner_reps = 10,
                     measure = "all", patch_size = 1,
                     quantiles = c(0.025, 0.5, 0.975),
                     silent = FALSE, pb_markdown = FALSE) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  if (is.null(loci)) {
    loci <- unique(project$data$raw$genetic_data$locus)
  }
  assert_vector_pos_int(loci, zero_allowed = FALSE)
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_single_pos_int(inner_reps, zero_allowed = FALSE)
  assert_single_string(measure)
  assert_in(measure, c("max_abs_grad", "mean_abs_grad", "variance", "all"))
  assert_single_pos_int(patch_size, zero_allowed = FALSE)
  assert_leq(patch_size, 5)
  if (!is.null(quantiles)) {
    assert_vector_bounded(quantiles)
  }
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # get which measures are turned on
  if (measure == "all") {
    measure <- c("max_abs_grad", "mean_abs_grad", "variance")
  }
  max_abs_grad_on <- ("max_abs_grad" %in% measure)
  mean_abs_grad_on <- ("mean_abs_grad" %in% measure)
  variance_on <- ("variance" %in% measure)
  
  # get adjacency list, taking patch size into account
  if (patch_size == 1) {
    adj_list <- project$maps$grid$adjacency
  } else {
    adj_mat <- project$maps$grid$adjacency %>% as.matrix()
    diag(adj_mat) <- 1
    for (i in 2:patch_size) {
      adj_mat <- adj_mat %*% adj_mat
    }
    adj_list <- apply(adj_mat, 1, function(x) which(x != 0))
  }
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$nu, .data$lambda, .data$omega, .data$gamma) %>%
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
  
  # get basic dimensions
  n_loci <- length(loci)
  n_alleles <- mapply(ncol, data_list)
  n_cells <- nrow(project$maps$grid$centroids)
  total_reps <- reps * inner_reps
  
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
  pb_dummy <- txtProgressBar(0, 1, initial = NA, style = 3)
  pb_post <- txtProgressBar(0, n_loci, initial = NA, style = 3)
  pb_null <- txtProgressBar(0, n_loci, initial = NA, style = 3)
  args_progress <- list(pb = pb_dummy,
                        pb_post = pb_post,
                        pb_null = pb_null)
  
  # create misc list
  args_misc <- list(silent = TRUE,
                    pb_markdown = TRUE)
  
  if (!silent) {
    message("Drawing from posterior distribution for each locus")
    update_progress(args_progress, "pb_post", 0, n_loci)
  }
  
  # get posterior surface
  post_max_abs_grad <- post_mean_abs_grad <- post_variance <- matrix(0, n_cells, total_reps)
  for (locus_i in seq_along(loci)) {
    
    # draw from predictive distribution via efficient C++ function
    sim_array <- predict_map_cpp(data_list[[locus_i]], mcmc_sample, dist_11, dist_12,
                                 dist_22, project$model$parameters, inner_reps,
                                 args_progress, args_functions, args_misc) %>%
      abind::abind(along = 2)
    
    #return(sim_array)
    
    # add measures to wombling posterior
    for (i in 1:n_cells) {
      adjacents <- adj_list[[i]]
      for (j in 1:n_alleles[locus_i]) {
        if (max_abs_grad_on || mean_abs_grad_on) {
          grad <- sweep(sim_array[adjacents,,j,drop = FALSE], 2, sim_array[i,,j], "-") %>% abs()
          if (max_abs_grad_on) {
            max_grad <- apply(grad, 2, max)
            post_max_abs_grad[i,] <- post_max_abs_grad[i,] + max_grad
          }
          if (mean_abs_grad_on) {
            mean_grad <- apply(grad, 2, mean)
            post_mean_abs_grad[i,] <- post_mean_abs_grad[i,] + mean_grad
          }
        }
        if (variance_on) {
          post_variance[i,] <- post_variance[i,] + apply(sim_array[adjacents,,j,drop = FALSE], 2, var)
        }
      }
    }
    
    if (!silent) {
      update_progress(args_progress, "pb_post", locus_i, n_loci)
    }
    
  }  # end loop over loci
  
  if (!silent) {
    message("Drawing from null distribution for each locus")
    update_progress(args_progress, "pb_null", 0, n_loci)
  }
  
  # get null surface
  null_max_abs_grad <- null_mean_abs_grad <- null_variance <- matrix(0, n_cells, total_reps)
  for (locus_i in seq_along(loci)) {
    
    # draw from null distribution via efficient C++ function
    sim_array <- null_map_cpp(data_list[[locus_i]], mcmc_sample, dist_11, dist_22,
                               project$model$parameters, inner_reps,
                               args_progress, args_functions, args_misc) %>%
      abind::abind(along = 2)
    
    #return(sim_array)
    
    # add measures to wombling null
    for (i in 1:n_cells) {
      adjacents <- adj_list[[i]]
      for (j in 1:n_alleles[locus_i]) {
        if (max_abs_grad_on || mean_abs_grad_on) {
          grad <- sweep(sim_array[adjacents,,j,drop = FALSE], 2, sim_array[i,,j], "-") %>% abs()
          if (max_abs_grad_on) {
            max_grad <- apply(grad, 2, max)
            null_max_abs_grad[i,] <- null_max_abs_grad[i,] + max_grad
          }
          if (mean_abs_grad_on) {
            mean_grad <- apply(grad, 2, mean)
            null_mean_abs_grad[i,] <- null_mean_abs_grad[i,] + mean_grad
          }
        }
        if (variance_on) {
          null_variance[i,] <- null_variance[i,] + apply(sim_array[adjacents,,j,drop = FALSE], 2, var)
        }
      }
    }
    
    if (!silent) {
      update_progress(args_progress, "pb_null", locus_i, n_loci)
    }
    
  }  # end loop over loci
  
  # get summaries and add to project
  if (max_abs_grad_on) {
    project$Wombling$max_abs_grad <- list(systemic_map = get_summaries(post_max_abs_grad / sum(n_alleles), quantiles = quantiles),
                                          null_map = get_summaries(null_max_abs_grad / sum(n_alleles), quantiles = quantiles),
                                          percentile_rank = rowMeans(post_max_abs_grad > null_max_abs_grad) * 100)
  }
  if (mean_abs_grad_on) {
    project$Wombling$mean_abs_grad <- list(systemic_map = get_summaries(post_mean_abs_grad / sum(n_alleles), quantiles = quantiles),
                                           null_map = get_summaries(null_mean_abs_grad / sum(n_alleles), quantiles = quantiles),
                                           percentile_rank = rowMeans(post_mean_abs_grad > null_mean_abs_grad) * 100)
  }
  if (variance_on) {
    project$Wombling$variance <- list(systemic_map = get_summaries(post_variance / sum(n_alleles), quantiles = quantiles),
                                      null_map = get_summaries(null_variance / sum(n_alleles), quantiles = quantiles),
                                      percentile_rank = rowMeans(post_variance > null_variance) * 100)
  }
  
  
  return(project)
}

#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'
#' @param project TODO
#' @param measure TODO
#' @param test_tail TODO
#' @param FDR TODO
#' @param silent TODO
#'
#' @export

Wombling_get_significant <- function(project, measure = "all", test_tail = "both",
                                     FDR = 0.05, silent = FALSE) {
  
  # check project
  assert_class(project, "genescaper_project")
  assert_single_string(measure)
  assert_in(measure, c("max_abs_grad", "mean_abs_grad", "variance", "all"))
  assert_single_string(test_tail)
  assert_in(test_tail, c("left", "right", "both"))
  assert_single_bounded(FDR)
  assert_single_logical(silent)
  
  if (!silent) {
    message("Calculating significance")
  }
  
  # loop through measures
  if (measure == "all") {
    measure <- c("max_abs_grad", "mean_abs_grad", "variance")
  }
  for (i in seq_along(measure)) {
    
    # get results into dataframe
    df_res <- data.frame(cell = seq_along(project$maps$grid$polygons),
                         prop = project$Wombling[[measure[i]]]$percentile_rank / 100)
    
    # calculate p-values
    if (test_tail == "left") {
      df_res$p <- df_res$prop
    } else if (test_tail == "right") {
      df_res$p <- 1.0 - df_res$prop
    } else if (test_tail == "both") {
      df_res$p <- ifelse(df_res$prop < 0.5, 2*df_res$prop, 2*(1 - df_res$prop))
    }
    
    # Bejamini and Yekutieli (2001) method for identifying significant results
    # while fixing the false descovery rate
    df_res$direction <- df_res$prop - 0.5
    which_signif <- Bejamini_Yekutieli(df_res, FDR)
    
    # add to project
    project$Wombling[[measure[i]]]$significance <- which_signif
  }
  
  return(project)
}



#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'
#' @param project TODO
#' @param loci TODO
#' @param reps TODO
#' @param quantiles TODO
#' @param silent TODO
#'
#' @export

get_mu_sigsq_credible <- function(project, loci, reps = 10,
                                  quantiles = c(0.025, 0.5, 0.975),
                                  silent = FALSE, pb_markdown = FALSE,
                                  true_sigsq) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_vector_pos_int(loci, zero_allowed = FALSE)
  assert_single_pos_int(reps, zero_allowed = FALSE)
  assert_vector_bounded(quantiles, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_logical(silent)
  assert_single_logical(pb_markdown)
  
  # sample parameters from posterior
  mcmc_sample <- project$model$MCMC$output %>%
    dplyr::filter(.data$phase == "sampling") %>%
    dplyr::select(.data$nu, .data$lambda, .data$omega, .data$gamma) %>%
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
  
  # TODO - remove
  allele_vec <- rep(loci, times = mapply(ncol, data_list))
  true_sigsq_list <- split(true_sigsq, f = allele_vec)
  
  # get distance between sampling sites
  dist_11 <- as.matrix(project$data$pairwise_measures$distance)
  
  # loop through loci
  mu_list <- sigsq_list <- list()
  for (i in seq_along(loci)) {
    
    # draw from predictive distribution via efficient C++ function
    z <- post_sigsq_mu(data_list[[i]], mcmc_sample, dist_11, true_sigsq_list[[i]])
    
    # get quantiles over mu
    mu_list[[i]] <- mapply(quantile, z$mu, list(probs = quantiles)) %>% t() %>%
      as.data.frame() %>%
      dplyr::mutate(locus = i,
                    allele = 1:ncol(data_list[[i]]),
                    .before = 1)
    
    # get quantiles over sigsq
    sigsq_list[[i]] <- mapply(quantile, z$sigsq, list(probs = quantiles)) %>% t() %>%
      as.data.frame() %>%
      dplyr::mutate(locus = i,
                    allele = 1:ncol(data_list[[i]]),
                    .before = 1)
    
  }
  
  # finalise outputs
  mu_df <-  dplyr::bind_rows(mu_list)
  sigsq_df <- dplyr::bind_rows(sigsq_list)
  
  return(list(mu = mu_df,
              sigsq = sigsq_df))
}
