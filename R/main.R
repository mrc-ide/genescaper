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
#' @param alpha_0 TODO
#' @param beta_0 TODO
#' @param phi_0 TODO
#' @param gamma_0 TODO
#'
#' @export

define_model <- function(project, alpha_0 = 1.0, beta_0 = 1.0, phi_0 = 0.0, gamma_0 = 0.1) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos(alpha_0, zero_allowed = FALSE)
  assert_single_pos(beta_0, zero_allowed = FALSE)
  assert_single_numeric(phi_0)
  assert_single_pos(gamma_0, zero_allowed = FALSE)
  
  # store within project
  project$model$parameters <- list(alpha_0 = alpha_0,
                                   beta_0 = beta_0,
                                   phi_0 = phi_0,
                                   gamma_0 = gamma_0)
  
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
    dplyr::select(.data$locus, .data$allele, .data$z) %>%
    dplyr::group_by(.data$locus, .data$allele) %>%
    dplyr::group_split() %>%
    lapply(function (x) x$z)
  
  # define parameters dataframe
  df_params <- drjacoby::define_params(name = "log_lambda", min = -Inf, max = Inf,
                                       name = "u", min = 0, max = 1)
  
  # source C++ likelihood and prior functions
  #Rcpp::sourceCpp("ignore/Cpp scripts/gp_model2.cpp")
  
  return(project)
}
