
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


#------------------------------------------------
#' @title TODO
#'
#' @description TODO
#'   
#' @param project a genescaper project, as produced by the
#'   \code{genescaper_project()} function. The project must have a grid already
#'   defined (see \code{?create_hex_grid}) as this grid will be used in
#'   simulation.
#' @param loci number of loci (assumed independent).
#' @param alleles number of alleles, assumed to be the same over all loci.
#' @param nu,lambda,omega the values of the model parameters used in simulation.
#'   These parameters together specify the spatial autocorrelation function.
#' @param sigsq_shape,sigsq_rate,mu_mean,gamma parameters of the
#'   normal-inverse-gamma prior on the mean and variance of transformed allele
#'   frequencies.
#'
#' @importFrom stats rnorm
#'
#' @export

sim_GRF <- function(project, loci, alleles, nu, lambda, omega,
                    sigsq_shape, sigsq_rate, mu_mean, gamma) {
  
  # check inputs
  assert_class(project, "genescaper_project")
  assert_single_pos_int(loci, zero_allowed = FALSE)
  assert_single_pos_int(alleles, zero_allowed = FALSE)
  assert_single_bounded(nu)
  assert_single_pos(lambda, zero_allowed = FALSE)
  assert_single_bounded(omega, left = 1, right = 2)
  assert_single_pos(sigsq_shape, zero_allowed = FALSE)
  assert_single_pos(sigsq_rate, zero_allowed = FALSE)
  assert_single_numeric(mu_mean)
  assert_single_pos(gamma, zero_allowed = FALSE)
  
  # define correlation matrix
  centroids <- project$maps$grid$centroids
  n_site <- nrow(centroids)
  d <- get_GC_distance(centroids$longitude, centroids$latitude)
  K <- diag(1 - nu, n_site) + nu * exp(-(d / lambda)^omega)
  
  # draw means and variances
  n_sim <- loci * (alleles - 1)
  sigsq <- 1 / rgamma(n_sim, shape = sigsq_shape, rate = sigsq_rate)
  mu <- rnorm(n_sim, mean = mu_mean, sd = sqrt(sigsq / gamma))
  df_params <- data.frame(locus = rep(1:loci, each = alleles - 1),
                          allele = 1:(alleles - 1),
                          mu = mu,
                          sigsq = sigsq)
  
  # simulate z-values all in one go
  z <- mvtnorm::rmvnorm(n_sim, sigma = K) %>%
    sweep(1, mu, "+") %>%
    sweep(1, sqrt(sigsq), "*")
  
  # apply transformation to each locus
  z_list <- split(as.data.frame(z), f = df_params$locus)
  sim_wide <- mapply(function(i) {
    ret <- transform_z_to_p(t(z_list[[i]]))
    colnames(ret) <- sprintf("allele_%s", seq_len(ncol(ret)))
    ret %>%
      as.data.frame() %>%
      dplyr::mutate(site_ID = seq_len(nrow(ret)),
                    locus = i)
  }, seq_along(z_list), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
  
  # get into long format
  sim_long <- sim_wide %>%
    tidyr::pivot_longer(-c(.data$site_ID, .data$locus), names_to = "allele", values_to = "freq") %>%
    dplyr::mutate(allele = as.numeric(as.factor(.data$allele)))
  
  return(list(data = sim_long,
              params = df_params))
}
