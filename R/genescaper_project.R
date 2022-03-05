
#------------------------------------------------
#' @title Define new GeneScapeR project
#'
#' @description TODO
#'
#' @export

genescaper_project <- function() {
  
  # create empty project
  project <- list(data = NULL,
                  model = NULL,
                  maps = NULL,
                  pairwise_predictions = NULL,
                  distance_predictions = NULL,
                  GeoMAPI = NULL)
  
  class(project) <- "genescaper_project"
  
  # return
  invisible(project)
}

#------------------------------------------------
# overload print() function
#' @method print genescaper_project
#' @export
print.genescaper_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function
#' @method summary genescaper_project
#' @export
summary.genescaper_project <- function(object, ...) {
  p <- object
  
  # if empty project
  if (all(mapply(is.null, p))) {
    message("(empty project)")
    invisible(object)
  }
  
  # print data properties
  if (!is.null(p$data)) {
    message("Data: (todo)")
    #n_demes <- length(p$epi_model_parameters$H)
    #message(sprintf("  demes: %s", n_demes))
    #message(sprintf("  H:\t %s", paste(p$epi_model_parameters$H, collapse = ", ")))
    #message(sprintf("  M:\t %s", paste(p$epi_model_parameters$M, collapse = ", ")))
    #message(sprintf("  seed infections: %s", paste(p$epi_model_parameters$seed_infections, collapse = ", ")))
    
    message("")
  }
  
  # print model parameters
  if (!is.null(p$model)) {
    message("Model: (todo)")
    
    message("")
  }
  
  invisible(object)
}

