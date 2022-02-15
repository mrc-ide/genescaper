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
