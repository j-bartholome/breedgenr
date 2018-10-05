# create_population -------------------------------------------------------
#'  Create or simulate a population based genotypic and phenotypic data
#'
#' \code{create_population} ceate a population object based on genotypic and phenotypic data
#'
#'#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param haplo an array of dimension p, n and 2 or a matrix of dimension p, 2n (NAs not  allowed for NW method).
#' @param pheno a data frame or a matrix of dimension n and y.
#'  The row names should identical to the column names of \code{haplo}.
#' @param map a data frame or a matrix of dimension p and 2.
#'  The row names should identical to the row names of \code{haplo}.
#' @param RR a numeric value corresponding de the recombination rate.
#' @param h2 a numeric value [0,1] or vector of length y.
#' @param simul logical scalar: TRUE (default) or FALSE.
#' @param method a character string indicating which method to use
#'   for the calibration of the genotype-phenotype relationship:
#'   "RF" (default), "SVR", or "NW".
#' @param nb_gen an integer indicating the number of generation
#'   to simulate from the current population. Default = 15
#' @param nb_ind_gen an integer or a vector of the same length \code{nb_gen}
#'   indicating the number of individual per generation. Default = 2000
#' @param nb_ind_P0 an integer indicating the final number of individuals in the population
#' @param folds an integer for the number of folds. Default  = 10.
#' @param n_rep an integer for the number of replicates. Default  = 10.
#' @param nb_cores the number of cores to use for parallel execution. Default = 1.
#' @param verbose logical scalar: TRUE (default) or FALSE
#' @param ... further arguments to be passed to or from methods.
#'
#' @return
#' A list :
#'  \itemize{
#'  \item{'haplo'}{ an array with the haplotype of the new population}
#'  \item{'geno'}{ a matrix with the genotype of the new population}
#' }
#' @export
#'
#' @seealso \code{\link{simul_geno}}, \code{\link{simul_P0}}, \code{\link{rf}}, \code{\link{svr}}
#'
#' @examples
#'
#' pop <- create_population(
#' haplo = rice1_x,
#' pheno = rice1_y,
#' map = rice1_map,
#' RR = 244000,
#' h2 =rep(0.3, 3),
#' nb_gen = 5,
#' nb_ind_gen = 200,
#' nb_ind_P0 = 100,
#' mtry = 1000
#' )
#'

create_population <-
  function(haplo = NULL,
           pheno = NULL,
           map = NULL,
           RR = NULL,
           h2 = NULL,
           simul = TRUE,
           method = c('RF', 'SVR', 'NW'),
           nb_gen = 15,
           nb_ind_gen = 2000,
           nb_ind_P0 = 500,
           folds = 10,
           n_rep = 10,
           nb_cores = 1,
           verbose = T,
           ...) {
    ellipsis <- list(...)
    stopifnot(is.logical(verbose), length(verbose) == 1)
    stopifnot(is.logical(simul), length(simul) == 1)
    if (is.null(haplo)) {
      stop("NULL object for 'haplo' and 'pheno' not allowed")
    }
    if (is.null(pheno)) {
      stop("NULL object for 'pheno' and 'pheno' not allowed")
    }
    if (is.null(map)) {
      stop("NULL object for 'map' not allowed")
    }
    if (!is.array(haplo) & !is.matrix(haplo)) {
      stop("'haplo' must be an array of dim p, n and 2, or
           a matrix of dim p, 2n not a ",
           class(haplo))
    } else {
      if (is.null(row.names(haplo)) | is.null(colnames(haplo))) {
        stop("'haplo' must have row (marker ID) and column (genotype ID) names")
      }
      if (!is.numeric(haplo)) {
        stop("'haplo' must be numeric")
      }
      if (length(which(is.na(haplo))) != 0) {
        stop("'haplo' must not contain missing values")
      }
      if (is.matrix(haplo)) {
        stopifnot(nrow(haplo) %% 2 == 0)
        haplo <-
          sapply(list(haplo[, seq(1, ncol(haplo), by = 2)],
                      haplo[, seq(2, ncol(haplo), by = 2)]),
                 identity, simplify = "array")
      }
    }
    pheno <- as.matrix(pheno)
    if (!is.matrix(pheno)) {
      stop("'pheno' must be a data.frame or a matrix not a ", class(pheno))
    } else {
      if (is.null(row.names(pheno))) {
        stop("'pheno' must have row names (genotype ID)")
      }
      if (!is.numeric(pheno)) {
        stop("'pheno' must be numeric")
      }
      if (length(which(is.na(pheno))) != 0) {
        stop("'pheno' must not contain missing values")
      }
    }
    if (!identical(sort(row.names(pheno)) , sort(colnames(haplo)))) {
      stop("genotype ID must be the same in 'haplo' and 'pheno'")
    } else{
      pheno <- pheno[colnames(haplo), , drop = F]
    }
    if (!is.data.frame(map))
      stop("'map' must be a data.frame, not a", class(map))
    if (!all(c("chromosome", "position") %in% colnames(map))) {
      stop("'map' must include the columns 'chromosome' and 'position' ")
    }
    stopifnot(is.numeric(map$position),
              identical(sort(rownames(map)) , sort(rownames(haplo))))
    map <- map[rownames(haplo),]
    if (simul == T) {
      stopifnot(is.numeric(h2))
      if (!length(h2) == ncol(pheno)) {
        stop(
          "'h2' must be a numeric vector [0,1] of the same length
          than the number of phenotype in 'pheno'"
        )
      }
      stopifnot(is.numeric(RR), length(RR) == 1)
      stopifnot(is.numeric(nb_gen), length(nb_gen) == 1)
      stopifnot(is.numeric(nb_ind_P0), length(nb_ind_P0) == 1)
      stopifnot(is.numeric(nb_cores), length(nb_cores) == 1)
      stopifnot(is.numeric(nb_ind_gen))
      if (!length(nb_ind_gen) == 1 &
          !length(nb_ind_gen) == nb_gen) {
        stop("'nb_ind_gen' must be either an integer or
             a vector of the same length than 'nb_gen'")
      }
      stopifnot(is.numeric(folds), length(folds) == 1)
      stopifnot(is.numeric(n_rep), length(n_rep) == 1)
      stopifnot(is.logical(verbose))
      method <- match.arg(method)
      if (method == 'RF') {
        if (is.null(ellipsis$mtry)) {
          stop("'mtry' argument is missing for method 'RF'")
        } else if (!is.numeric(ellipsis$mtry)) {
          stop("'mtry' argument for method 'RF' must be an integer")
        }
      } else if (method == "SVR") {
        if (is.null(ellipsis$C)) {
          stop("'C' argument is missing for method 'SVR'")
        } else if (!is.numeric(ellipsis$C)) {
          stop("'C' argument for method 'SVR' must be an integer")
        }
        if (!is.numeric(ellipsis$sigma)) {
          stop("'sigma' argument for method 'SVR' must be numeric")
        }
      } else if (method == "NW") {
        if (is.null(ellipsis$h)) {
          # voir comment faire pour h utiliser sigest?
        } else if (!is.numeric(ellipsis$h)) {
          stop("'h' argument for method 'NW' must be an integer")
        }
        kernel_type <-
          match.arg(ellipsis$kernel_type, c('rho', 'gaussian', 'aggregate'))
      }
      pop <- simul_P0(
        haplo = haplo,
        pheno = pheno,
        map = map,
        RR = RR,
        h2 = h2,
        nb_gen = nb_gen,
        nb_ind_gen = nb_ind_gen,
        nb_ind_P0 = nb_ind_P0,
        folds = folds,
        n_rep = n_rep,
        method = method,
        nb_cores = nb_cores,
        verbose = verbose,
        ...
      )
      } else {
        geno <- haplo[, , 1] + haplo[, , 2]
        geno <- t(geno)
        pop <- list(
          haplo = haplo,
          geno = geno,
          phen_val = setNames(split(as.matrix(pheno), col(
            pheno)), colnames(pheno)),
          map = map
        )
      }
    return(pop)
    }

