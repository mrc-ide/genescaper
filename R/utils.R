#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
# density of logit-normal distribution
#' @importFrom stats dnorm
dlogitnorm <- function(p, raw_mu = 0.0, raw_sigsq = 1.0, return_log = TRUE) {
  ret <- dnorm(log(p) - log(1.0 - p), mean = raw_mu, sd = sqrt(raw_sigsq), log = TRUE) - log(p) - log(1.0 - p)
  if (!return_log) {
    ret <- exp(ret)
  }
  return(ret)
}

#------------------------------------------------
#' @title Calculate great circle distance between pairwise coordinates
#'
#' @description Get great circle distance between spatial points defined by a
#'   vector of longitudes and latitudes. Distances are returned in a pairwise
#'   distance matrix.
#' 
#' @param long,lat vector of longitudes and latitudes.
#'
#' @importFrom  stats as.dist
#' @export

get_GC_distance <- function(long, lat) {
  
  # check inputs
  assert_vector_numeric(long)
  assert_vector_numeric(lat)
  assert_same_length(long, lat)
  assert_bounded(long, left = -180, right = 180)
  assert_bounded(lat, left = -90, right = 90)
  
  # calculate distance matrix
  ret <- apply(cbind(long, lat), 1, function(y) {lonlat_to_bearing(long, lat, y[1], y[2])$gc_dist})
  ret <- as.dist(ret, upper = FALSE)
  
  return(ret)
}

#------------------------------------------------
#' @title Calculate great circle distance and bearing between coordinates
#'
#' @description Calculate great circle distance and bearing between spatial
#'   coordinates, defined by longitude and latitude of both origin and
#'   destination points.
#'
#' @param origin_lon,origin_lat the origin longitude and latitude.
#' @param dest_lon,dest_lat the destination longitude and latitude
#'
#' @export
#' @examples
#' # one degree longitude should equal approximately 111km at the equator
#' lonlat_to_bearing(0, 0, 1, 0)

lonlat_to_bearing <- function(origin_lon, origin_lat, dest_lon, dest_lat) {
  
  # check inputs
  assert_vector_numeric(origin_lon)
  assert_vector_numeric(origin_lat)
  assert_vector_numeric(dest_lon)
  assert_vector_numeric(dest_lat)
  
  # convert input arguments to radians
  origin_lon <- origin_lon * 2 * pi / 360
  origin_lat <- origin_lat * 2 * pi / 360
  dest_lon <- dest_lon * 2 * pi / 360
  dest_lat <- dest_lat * 2 * pi / 360
  
  # get change in lon
  delta_lon <- dest_lon - origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat),
                   cos(origin_lat)*sin(dest_lat) - sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or 
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  
  # convert bearing from radians to degrees measured clockwise from due north,
  # and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing * 360 / (2 * pi)
  bearing <- (bearing + 360) %% 360
  earth_rad <- 6371
  gc_dist <- earth_rad*gc_angle
  
  # return list
  ret <-list(bearing = bearing,
             gc_dist = gc_dist)
  return(ret)
}

#------------------------------------------------
# Transform frequencies (in interval [0,1]) to continuous z values (in interval
# [-Inf, Inf])
#' @noRd
transform_p_to_z <- function(genetic_data) {
  
  # check inputs
  assert_dataframe(genetic_data)
  assert_in(c("site_ID", "locus", "allele", "freq"), names(genetic_data))
  
  # get adjusted allele frequencies by applying stick breaking correction. For
  # each allele, p_stick is calculated by dividing the allele frequency by the
  # amount of unit interval that remains once previous frequencies have been
  # taken into account. This creates (J - 1) ostensibly independent frequencies,
  # where J is the number of alleles. The final p_stick always equals 1 and so
  # this allele can be dropped. Iin other words, there are J - 1 degrees of
  # freedom and so we are reducing from J dependent frequencies to J-1
  # independent frequencies. z values are then calculated as logit(p_stick).
  ret <- genetic_data %>%
    dplyr::group_by(.data$site_ID, .data$locus) %>%
    dplyr::summarise(J = length(.data$allele),
                     allele = .data$allele[-.data$J],
                     stick_remaining = 1 - cumsum(.data$freq[-.data$J]) + .data$freq[-.data$J],
                     p_stick = .data$freq[-.data$J] / .data$stick_remaining,
                     z = log(.data$p_stick) - log(1 - .data$p_stick)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$J, -.data$stick_remaining, -.data$p_stick)
  
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}
