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
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x - x_max)))
  return(ret)
}

#------------------------------------------------
# improved version of log_sum() that gives higher accuracy (less sensitive to
# underflow) in situations when the largest value of x is 0, and the second
# largest value is much more negative. For example, log_sum(c(0, -40)) returns
# 0, whereas log_sum2(c(0, -40)) returns a small positive value.
#' @noRd
log_sum2 <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  w <- which.max(x)
  if ((x[w] == 0) && (max(x[-w]) < -36)) {
    ret <- sum(exp(x[-w] - x[w]))
  } else {
    ret <- log_sum(x)
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
#' @param lon,lat vector of longitudes and latitudes.
#'
#' @importFrom  stats as.dist
#' @export

get_GC_distance <- function(lon, lat) {
  
  # check inputs
  assert_vector_numeric(lon)
  assert_vector_numeric(lat)
  assert_same_length(lon, lat)
  assert_bounded(lon, left = -180, right = 180)
  assert_bounded(lat, left = -90, right = 90)
  
  # calculate distance matrix
  ret <- apply(cbind(lon, lat), 1, function(y) {lonlat_to_bearing(lon, lat, y[1], y[2])$gc_dist})
  
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
  # this allele can be dropped. In other words, there are J - 1 degrees of
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
# The reverse transform of transform_p_to_z(). Takes a matrix of z-values with
# sites in rows and alleles in columns.
#' @noRd
transform_z_to_p <- function(z) {
  
  # check inputs
  assert_matrix_numeric(z)
  
  # get basic dimensions
  n_alleles <- ncol(z)
  
  # logistic transform
  q <- 1 / (1 + exp(-z))
  
  # apply stick-breaking transform
  p <- matrix(NA, nrow(q), n_alleles + 1)
  p[,1] <- q[,1]
  stick_remaining <- 1 - p[,1]
  if (n_alleles > 1) {
    for (i in 2:n_alleles) {
      p[,i] <- stick_remaining * q[,i]
      stick_remaining <- stick_remaining - p[,i]
    }
  }
  p[,n_alleles + 1] <- stick_remaining
  
  return(p)
}

#------------------------------------------------
# Equivalent to transform_z_to_p(), but evaluated in log space (returning
# log(p)), and operates on a 3D array with the following dimensions: 1) demes,
# 2) reps, 3) alleles. Uses some approximations for extreme large or small z to
# avoid underflow issues that otherwise could result in logistic transform
# returning 0 or 1, which should not be possible for finite z.
#' @noRd
transform_z_to_logp <- function(z) {
  
  # check inputs
  assert_length(dim(z), 3)
  
  # get basic dimensions
  n_demes <- dim(z)[1]
  n_reps <- dim(z)[2]
  n_alleles <- dim(z)[3]
  
  # log-logistic transform
  log_q <- -log(1 + exp(-z))
  
  # use approximations for extreme values to avoid underflow
  log_q[z < -100] <- z[z < -100]
  log_q[z > 100] <- -exp(-z[z > 100])
  
  # apply stick-breaking transform
  log_p <- log_q
  for (i in 2:n_alleles) {
    if (i == 2) {
      log_sum_p <- log_p[,,1]
    } else {
      log_sum_p <- apply(log_p[,,1:(i - 1), drop = FALSE], c(1, 2), log_sum2)
    }
    log_stick_remaining <- log(1 - exp(log_sum_p))
    w <- which(log_sum_p < -36)
    log_stick_remaining[w] <- -exp(log_sum_p[w])
    w <- which((log_sum_p < 0) & (log_sum_p > -1e-15))
    log_stick_remaining[w] <- log(-log_sum_p[1])
    if (i == n_alleles) {
      log_p[,,i] <- log_stick_remaining
    } else {
      log_p[,,i] <- log_q[,,i] + log_stick_remaining
    }
  }
  
  return(log_p)
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

#------------------------------------------------
# if a single value is provided then expand to a vector of length n
#' @noRd
force_vector <- function(x, n) {
  if (length(x) == 1) {
    return(rep(x,n))
  } else {
    return(x)
  }
}

#------------------------------------------------
#' @title Calculate ellipse polygon coordinates from foci and eccentricity
#'
#' @description TODO calculate ellipse polygon coordinates from foci and eccentricity.
#'
#' @param f1 x- and y-coordinates of the first focus.
#' @param f2 x- and y-coordinates of the second focus.
#' @param ecc eccentricity of the ellipse, defined as half the distance between
#'   foci divided by the semi-major axis. We can say \eqn{e = sqrt{1 -
#'   b^2/a^2}}, where \eqn{e} is the eccentricity, \eqn{a} is the length of the
#'   semi-major axis, and \eqn{b} is the length of the semi-minor axis.
#'   Eccentricity ranges from 0 (perfect circle) to 1 (straight line between
#'   foci), although eccentricity of 0 is not allowed in this function as it
#'   would lead to an infinitely large circle.
#' @param n number of points in polygon.
#'
#' @export

get_ellipse <- function(f1 = c(-3, 0), f2 = c(3, 0), ecc = 0.8, n = 100) {
  
  # check inputs
  assert_vector_numeric(f1)
  assert_length(f1, 2)
  assert_vector_numeric(f2)
  assert_length(f2, 2)
  assert_single_bounded(ecc, inclusive_left = FALSE)
  assert_bounded(ecc, inclusive_left = FALSE)
  assert_single_pos_int(n)
  
  # define half-distance between foci (c), semi-major axis (a) and semi-minor
  # axis(b)
  c <- 0.5 * sqrt(sum((f2 - f1)^2))
  a <- c / ecc
  b <- sqrt(a^2 - c^2)
  
  # define slope of ellipse (alpha) and angle of points from centre (theta)
  alpha <- atan2(f2[2] - f1[2], f2[1] - f1[1])
  theta <- seq(0, 2*pi, l = n + 1)
  
  # define x and y coordinates
  x <- (f1[1] + f2[1]) / 2 + a*cos(theta)*cos(alpha) - b*sin(theta)*sin(alpha)
  y <- (f1[2] + f2[2]) / 2 + a*cos(theta)*sin(alpha) + b*sin(theta)*cos(alpha)
  
  # ensure ellipse closes perfectly
  x[n + 1] <- x[1]
  y[n + 1] <- y[1]
  
  # return as dataframe
  return(data.frame(x = x, y = y))
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# calculate pairwise Gst from a matrix or 3D array of allele frequencies. The
# dimensions must be either 1) demes, 2) alleles for a matrix, or 1) demes, 2)
# reps, 3) alleles for an array.
#' @noRd
calc_pairwise_Gst <- function(x) {
  
  # split based on matrix vs array
  if (length(dim(x)) == 2) {
    
    # get basic dimensions
    demes <- nrow(x)
    alleles <- ncol(x)
    
    # get all pairs of demes
    deme_pairs <- cbind(rep(1:(demes - 1), times = (demes - 1):1),
                        unlist(mapply(function(i) i:demes, 2:demes)))
    
    # calculate Gst for all pairs
    x_bar <- 0.5 * (x[deme_pairs[,1],] + x[deme_pairs[,2],])
    J_t <- rowSums(x_bar^2)
    J_i <- rowSums(x^2)
    J_s <- 0.5 * (J_i[deme_pairs[,1]] + J_i[deme_pairs[,2]])
    
    # deal with underflow issues by restricting range of J_s and J_t. Note that
    # allele frequencies of exactly 0 or 1 should not be possible under the
    # model, therefore both J_s and J_t should be restricted to be <1. They
    # should also be restricted to be >0 because this would require the sum of
    # infinitely many values. Finally, it should not be possible for J_t to
    # exceed J_s.
    J_s[J_s < 1e-10] <- 1e-10
    J_s[J_s > 1.0 - 1e-10] <- 1 - 1e-10
    J_t[J_t < 1e-10] <- 1e-10
    J_t[J_t > 1.0 - 1e-10] <- 1 - 1e-10
    J_s[J_s < J_t] <- J_t[J_s < J_t]
    
    Gst <- (J_s - J_t) / (1 - J_t)
    
  } else if (length(dim(x)) == 3) {
    
    # get basic dimensions
    demes <- dim(x)[1]
    reps <- dim(x)[2]
    alleles <- dim(x)[3]
    
    # get all pairs of demes
    deme_pairs <- cbind(rep(1:(demes - 1), times = (demes - 1):1),
                        unlist(mapply(function(i) i:demes, 2:demes)))
    
    # calculate Gst for all pairs
    x_bar <- 0.5 * (x[deme_pairs[,1],,,drop = FALSE] + x[deme_pairs[,2],,,drop = FALSE])
    J_t <- J_i <- 0
    for (i in seq_len(alleles)) {
      J_t <- J_t + x_bar[,,i]^2
      J_i <- J_i + x[,,i]^2
    }
    if (reps == 1) {
      J_t <- matrix(J_t)
      J_i <- matrix(J_i)
    }
    J_s <- 0.5 * (J_i[deme_pairs[,1],,drop = FALSE] + J_i[deme_pairs[,2],,drop = FALSE])
    
    # deal with underflow issues by restricting range of J_s and J_t. Note that
    # allele frequencies of exactly 0 or 1 should not be possible under the
    # model, therefore both J_s and J_t should be restricted to be <1. They
    # should also be restricted to be >0 because this would require the sum of
    # infinitely many values. Finally, it should not be possible for J_t to
    # exceed J_s.
    J_s[J_s < 1e-10] <- 1e-10
    J_s[J_s > 1.0 - 1e-10] <- 1 - 1e-10
    J_t[J_t < 1e-10] <- 1e-10
    J_t[J_t > 1.0 - 1e-10] <- 1 - 1e-10
    J_s[J_s < J_t] <- J_t[J_s < J_t]
    
    Gst <- (J_s - J_t) / (1 - J_t)
    
  } else {
    stop("input must be matrix or array")
  }
  
  return(Gst)
}

# -----------------------------------
# calculate pairwise Jost's D from a matrix or 3D array of allele frequencies. The
# dimensions must be either 1) demes, 2) alleles for a matrix, or 1) demes, 2)
# reps, 3) alleles for an array.
#' @noRd
calc_pairwise_D <- function(x) {
  
  # split based on matrix vs array
  if (length(dim(x)) == 2) {
    
    # get basic dimensions
    demes <- nrow(x)
    alleles <- ncol(x)
    
    # get all pairs of demes
    deme_pairs <- cbind(rep(1:(demes - 1), times = (demes - 1):1),
                        unlist(mapply(function(i) i:demes, 2:demes)))
    
    # calculate D for all pairs
    x_bar <- 0.5 * (x[deme_pairs[,1],] + x[deme_pairs[,2],])
    J_t <- rowSums(x_bar^2)
    J_i <- rowSums(x^2)
    J_s <- 0.5 * (J_i[deme_pairs[,1]] + J_i[deme_pairs[,2]])
    
    # deal with underflow issues by restricting range of J_s and J_t. Note that
    # allele frequencies of exactly 0 or 1 should not be possible under the
    # model, therefore both J_s and J_t should be restricted to be <1. They
    # should also be restricted to be >0 because this would require the sum of
    # infinitely many values. Finally, it should not be possible for J_t to
    # exceed J_s.
    J_s[J_s < 1e-10] <- 1e-10
    J_s[J_s > 1.0 - 1e-10] <- 1 - 1e-10
    J_t[J_t < 1e-10] <- 1e-10
    J_t[J_t > 1.0 - 1e-10] <- 1 - 1e-10
    J_s[J_s < J_t] <- J_t[J_s < J_t]
    
    D <- (J_s - J_t) / J_s * demes / (demes - 1)
    
  } else if (length(dim(x)) == 3) {
    
    # get basic dimensions
    demes <- dim(x)[1]
    reps <- dim(x)[2]
    alleles <- dim(x)[3]
    
    # get all pairs of demes
    deme_pairs <- cbind(rep(1:(demes - 1), times = (demes - 1):1),
                        unlist(mapply(function(i) i:demes, 2:demes)))
    
    # calculate Gst for all pairs
    x_bar <- 0.5 * (x[deme_pairs[,1],,,drop = FALSE] + x[deme_pairs[,2],,,drop = FALSE])
    J_t <- J_i <- 0
    for (i in seq_len(alleles)) {
      J_t <- J_t + x_bar[,,i]^2
      J_i <- J_i + x[,,i]^2
    }
    if (reps == 1) {
      J_t <- matrix(J_t)
      J_i <- matrix(J_i)
    }
    J_s <- 0.5 * (J_i[deme_pairs[,1],,drop = FALSE] + J_i[deme_pairs[,2],,drop = FALSE])
    
    # deal with underflow issues by restricting range of J_s and J_t. Note that
    # allele frequencies of exactly 0 or 1 should not be possible under the
    # model, therefore both J_s and J_t should be restricted to be <1. They
    # should also be restricted to be >0 because this would require the sum of
    # infinitely many values. Finally, it should not be possible for J_t to
    # exceed J_s.
    J_s[J_s < 1e-10] <- 1e-10
    J_s[J_s > 1.0 - 1e-10] <- 1 - 1e-10
    J_t[J_t < 1e-10] <- 1e-10
    J_t[J_t > 1.0 - 1e-10] <- 1 - 1e-10
    J_s[J_s < J_t] <- J_t[J_s < J_t]
    
    D <- (J_s - J_t) / J_s * demes / (demes - 1)
    
  } else {
    stop("input must be matrix or array")
  }
  
  return(D)
}

# -----------------------------------
# takes a series of n-choose-2 pairwise values and aranges them in an n-by-n
# matrix, in which values are reflected on the diagonal
#' @noRd
pairwise_to_mat <- function(x) {
  
  # get required dimension of matrix
  n <- 0.5 + 0.5 * sqrt(1 + 8 * length(x))
  
  # arrange in matrix
  ret <- matrix(0, n, n)
  ret[lower.tri(ret)] <- x
  ret[upper.tri(ret)] <- t(ret)[upper.tri(ret)]
  
  ret
}

# -----------------------------------
# Bejamini and Yekutieli (2001) method for identifying significant results while
# fixing the false descovery rate. df_res must be a dataframe with columns:
# cell, p, direction.
#' @noRd
Bejamini_Yekutieli <- function(df_res, FDR) {
  
  # sort in order of increasing p
  df_res <- df_res[order(df_res$p),]
  
  # Bejamini and Yekutieli (2001) method for identifying significant results
  # while fixing the false descovery rate
  df_res$BY <- FDR * seq_along(df_res$p) / nrow(df_res)
  which_lower <- which_upper <- integer()
  if (any(df_res$p <= df_res$BY, na.rm = TRUE)) {
    
    w <- which(df_res$p <= df_res$BY)
    which_upper <- df_res$cell[w][df_res$direction[w] > 0]
    which_lower <- df_res$cell[w][df_res$direction[w] <= 0]
  }
  
  return(list(upper = which_upper,
              lower = which_lower))
}

# -----------------------------------
# get mean, sd and quantiles over rows of matrix
#' @noRd
get_summaries <- function(x, quantiles = NULL) {
  ret <- list(mean = rowMeans(x),
              sd = apply(x, 1, sd))
  if (!is.null(quantiles)) {
    quants <- list()
    for (i in seq_along(quantiles)) {
      quants[[i]] <- apply(x, 1, quantile, probs = quantiles[i])
    }
    names(quants) <- sprintf("%s%%", round(quantiles * 100, digits = 1))
    ret$quantiles <- quants
  }
  return(ret)
}

#------------------------------------------------
# pass in a series of polygons (class sfc_POLYGON). Expand by a buffer distance
# d and merge polygons
#' @noRd
get_merged_poly <- function(hex_polys, d = 0.1) {
  
  # expand polygons
  ret <- st_buffer(hex_polys, d)
  
  # merge polygons
  ret <- sf::st_union(sf::st_sf(ret))
  
  # undo expansion
  ret <- st_buffer(ret, -d)
  
  return(ret)
}
