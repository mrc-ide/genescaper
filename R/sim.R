
#------------------------------------------------
#' @title Get distance between points taking into account barriers
#'
#' @description Given a set of lon/lat coordinates and a list of barriers in the
#'   form of polygons, returns the "effective distance" between points, which is
#'   calculated as the great-circle distance with a penalty applied if the line
#'   intersects a barrier. The exact way in which barriers modify distances can
#'   be varied (see \code{barrier_method} argument).
#'
#' @param node_lon,node_lat longitudes and latitudes of nodes.
#' @param barrier_list list of polygons representing barriers. Each element of
#'   the list must be a dataframe with columns \code{longitude} and
#'   \code{latitude} specifying the coordinates of points that make up the
#'   polygon. Polygons must be complete rings, meaning the final row of the
#'   dataframe must equal the first row.
#' @param barrier_penalty penalty values of each barrier. If a single value is
#'   provided then this value will be used for all barriers.
#' @param barrier_method the method by which penalties are applied:
#'   \enumerate{
#'     \item{compare pairwise lines to barriers. If the line intersects then add
#'     a fixed \code{barrier_penalty} to the spatial distance.}
#'     \item{compare pairwise lines to barriers. Calculate the intersection of
#'     the two, multiply this by the \code{barrier_penalty} and add to the
#'     spatial distance. For example, a \code{barrier_penalty} of 1 would mean
#'     there is double the "friction" when moving through a barrier.}
#'     \item{compare pairwise ellipses to barriers. Calculate the intersection
#'     area of the two, multiply this by the \code{barrier_penalty} and add to
#'     the spatial distance.}
#'   }
#' @param max_barrier_range edges that are longer than this distance are
#'   unaffected by any barriers. Makes it possible to model barriers that only
#'   apply locally.
#' @param eccentricity eccentricity of ellipses (only used under
#'   \code{barrier_method = 3}).
#' @param noise_sd standard deviation of Gaussian noise added to all distances
#'   (after the application of barriers).
#' @param ellipse_points number of points that make up an ellipse (only used
#'   under \code{barrier_method = 3}).
#'
#' @import sf
#' @importFrom stats dist rnorm
#' @export

get_barrier_intersect <- function(node_lon,
                                  node_lat,
                                  barrier_list = list(),
                                  barrier_penalty = numeric(),
                                  barrier_method = 1,
                                  max_barrier_range = Inf,
                                  eccentricity = 0.9,
                                  noise_sd = 0,
                                  ellipse_points = 20) {
  
  # check inputs
  assert_vector_numeric(node_lon)
  assert_vector_numeric(node_lat)
  assert_same_length(node_lon, node_lat)
  assert_list(barrier_list)
  nb <- length(barrier_list)
  if (nb > 0) {
    for (i in 1:nb) {
      assert_dataframe(barrier_list[[i]])
      assert_in(c("longitude", "latitude"), names(barrier_list[[i]]))
      assert_eq(barrier_list[[i]][1,], barrier_list[[i]][nrow(barrier_list[[i]]),], 
                message = "barrier polygons must be closed, i.e. the last node coordinate equals the first")
    }
  }
  assert_vector_numeric(barrier_penalty)
  assert_single_pos_int(barrier_method)
  assert_in(barrier_method, 1:3)
  assert_single_pos(max_barrier_range, zero_allowed = TRUE)
  assert_single_bounded(eccentricity, inclusive_left = FALSE)
  assert_single_pos(noise_sd, zero_allowed = TRUE)
  assert_single_pos_int(ellipse_points, zero_allowed = FALSE)
  
  # force barrier_penalty to vector
  barrier_penalty <- force_vector(barrier_penalty, length(barrier_list))
  assert_same_length(barrier_penalty, barrier_list)
  
  # create mask for ignoring edges greater then
  distance_mask <- 1
  if (is.finite(max_barrier_range)) {
    d <- as.vector(get_GC_distance(node_lon, node_lat))
    distance_mask <- (d < max_barrier_range)
  }
  
  # apply barrier penalties
  intersect_penalty <- 0
  if ((nb > 0) & any(barrier_penalty != 0)) {
    
    # convert barrier list to st_polygon
    poly_list <- list()
    for (i in 1:length(barrier_list)) {
      poly_list[[i]] <- sf::st_polygon(list(as.matrix(barrier_list[[i]])))
    }
    
    # get node coordinates in matrix
    node_mat <- cbind(node_lon, node_lat)
    
    # if comparing lines
    if (barrier_method %in% c(1, 2)) {
      
      # create all pairwise sf_linestring between nodes
      line_list <- list()
      n_node <- length(node_lon)
      i2 <- 0
      for (i in 1:(n_node-1)) {
        for (j in (i + 1):n_node) {
          i2 <- i2 + 1
          line_list[[i2]] <- sf::st_linestring(node_mat[c(i,j),])
        }
      }
      
      # convert lines and polys to st_sfc
      line_sfc <- sf::st_sfc(line_list)
      poly_sfc <- sf::st_sfc(poly_list)
      
      # get boolean intersection matrix
      intersect_mat <- as.matrix(sf::st_intersects(line_sfc, poly_sfc))
      
      # convert to (great circle) length of intersection if using method 2
      if (barrier_method == 2) {
        intersect_mat[intersect_mat == TRUE] <- mapply(function(x) {
          
          # sum great circle distances of all intersection lenfths
          if ("LINESTRING" %in% class(x)) {
            get_GC_distance(as.matrix(x)[1,], as.matrix(x)[2,])[1]
          } else if ("MULTILINESTRING" %in% class(x)) {
            sum(mapply(function(y) {
              get_GC_distance(y[1,], y[2,])[1]
            }, x))
          } else {
            warning("cannot calculate distance: sfc type not recognised")
          }
          
        }, sf::st_intersection(line_sfc, poly_sfc))
      }
    }
    
    # mask out edges that are beyond limit distance
    intersect_mat <- sweep(intersect_mat, 1, distance_mask, "*")
    
    # if comparing ellipse
    if (barrier_method == 3) {
      
      # create all pairwise ellipses between nodes
      ell_list <- list()
      n_node <- length(node_lon)
      i2 <- 0
      for (i in 1:(n_node-1)) {
        for (j in (i+1):n_node) {
          i2 <- i2 + 1
          ell_df <- get_ellipse(f1 = node_mat[i,], f2 = node_mat[j,], ecc = eccentricity, n = ellipse_points)
          ell_list[[i2]] <- sf::st_polygon(list(as.matrix(ell_df)))
        }
      }
      
      # convert ellipses and polys to st_sfc
      ell_sfc <- sf::st_sfc(ell_list)
      poly_sfc <- sf::st_sfc(poly_list)
      
      # get boolean intersection matrix
      intersect_mat <- as.matrix(sf::st_intersects(ell_sfc, poly_sfc))
      
      # mask out ellipses that are beyond limit distance
      intersect_mat <- sweep(intersect_mat, 1, distance_mask, "*")
      
      # convert to area of intersection
      intersect_areas <- mapply(function(x) {
        sf::st_area(x)
      }, sf::st_intersection(ell_sfc, poly_sfc))
      intersect_mat[intersect_mat == TRUE] <- intersect_areas[intersect_mat == TRUE]
      
    }
    
    # apply penalty
    intersect_penalty <- rowSums(sweep(intersect_mat, 2, barrier_penalty, '*'))
    
  }  # end apply barrier penalties
  
  # get pairwise distance plus penalty
  d <- get_GC_distance(node_lon, node_lat) + intersect_penalty
  d <- d + rnorm(length(d), sd = noise_sd)
  
  # return matrix
  return(as.matrix(d))
}

#------------------------------------------------
#' @title Simulate allele frequencies from Wright-Fisher model
#'
#' @description Simulate Wright-Fisher evolution in a series of partially
#'   connected demes.
#'
#' @details Assumes a haploid population and independent loci (no linkage
#'   disequilibrium). Implements a finite-alleles mutation model with equal
#'   chance of mutating from any allele to any other. Migration is implemented
#'   by proposing random swaps of individuals between demes, thereby ensuring
#'   population sizes remain constant over time. For this reason, \code{N} must
#'   be the same for all demes.
#'   
#' @param N number of individuals per deme. Must be the same for all
#'   demes.
#' @param L number of loci (assumed independent).
#' @param alleles number of alleles. Can be a single number for all loci or a
#'   vector of length \code{L}.
#' @param mu mutation rate. Assumes finite-alleles model, with equal chance of
#'   mutating from any allele to any other.
#' @param mig_mat migration matrix specifying the per-generation probability of
#'   an individual migrating from any deme (in rows) to any other deme (in
#'   columns).
#' @param t_out vector of times at which results will be output.
#' @param initial_method,initial_params method of initialising allele
#'   frequencies, and parameters that are used in initialisation. There are two
#'   possible options:
#'   \enumerate{
#'   \item Each deme has allele frequencies drawn independently from a symmetric
#'   Dirichlet(theta/k) distribution, where \eqn{theta = 2*N*mu} and k is the
#'   number of alleles. This is the analytical equilibrium distribution under
#'   the model if there was no migration between demes.
#'     \item All demes have the same initial allele frequencies, which are drawn
#'     once from a Dirichlet(alpha_1, ..., alpha_k) distribution, where the
#'     alpha parameters are input as \code{initial_params}. This can be a vector
#'     if the same number of alleles is used over all loci, or a list of vectors
#'     over loci to accommodate varying numbers of alleles.
#'   }
#' @param silent if \code{TRUE} then suppress output to console.
#'
#' @importFrom utils txtProgressBar
#' @export

sim_wrightfisher <- function(N, L, alleles, mu, mig_mat, t_out,
                             initial_method = 1, initial_params = NULL, silent = FALSE) {
  
  # check inputs
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos_int(L, zero_allowed = FALSE)
  assert_vector_pos_int(alleles, zero_allowed = FALSE)
  assert_gr(alleles, 1)
  if (length(alleles) > 1) {
    assert_length(alleles, L)
  }
  assert_single_bounded(mu)
  assert_symmetric_matrix(mig_mat)
  assert_bounded(mig_mat)
  if (!all.equal(rowSums(mig_mat), rep(1, nrow(mig_mat)), check.attributes = FALSE)) {
    stop("every row of mig_mat must sum to 1")
  }
  assert_vector_pos_int(t_out, zero_allowed = TRUE)
  assert_in(initial_method, c(1, 2))
  if (initial_method == 2) {
    if (length(alleles) == 1) {
      assert_vector(initial_params)
      assert_length(initial_params, alleles)
    } else {
      assert_list(initial_params)
      assert_length(initial_params, L)
      for (i in 1:L) {
        assert_length(initial_params[[i]], alleles[i])
      }
    }
  }
  assert_single_logical(silent)
  
  # process some inputs
  if (length(alleles) == 1) {
    alleles <- rep(alleles, L)
  }
  if (initial_method == 2) {
    if (!is.list(initial_params)) {
      initial_params = replicate(L, initial_params, simplify = FALSE)
    }
  }
  
  # get number of demes from dimensions of migration matrix
  K <- ncol(mig_mat)
  
  # make argument list
  args <- list(N = N,
               K = K,
               L = L,
               alleles = alleles,
               mu = mu,
               mig_mat = matrix_to_rcpp(mig_mat),
               t_out = t_out,
               initial_method = initial_method,
               initial_params = initial_params,
               silent = silent)
  
  # create progress bars
  pb <- txtProgressBar(min = 0, max = max(c(1, t_out)), initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # functions to pass to C++
  args_functions <- list(update_progress = update_progress)
  
  # run efficient C++ function
  output_raw <- sim_wrightfisher_cpp(args, args_functions, args_progress)
  
  # process output
  output_processed <- mapply(function(t_i) {
    mapply(function(k) {
      data.frame(time = t_out[t_i],
                 deme = k,
                 locus = rep(seq_len(L), times = alleles),
                 allele = unlist(lapply(alleles, seq_len)),
                 count = unlist(output_raw$allele_counts[[t_i]][[k]]))
    }, seq_len(K), SIMPLIFY = FALSE)
  }, seq_along(t_out), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
  
  return(output_processed)
}
