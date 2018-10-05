# simul_P0 ----------------------------------------------------------------
#'  Simulate a population
#'
#' \code{simul_P0} simulate a population based on genotypic and phenotypic data
#'
#' @param haplo an array of dimension p, n and 2 or a matrix of dimension p, 2n (NAs not  allowed for NW method).
#' @param pheno a data frame or a matrix of dimension n and y.
#'  The row names should identical to the column names of \code{haplo}.
#' @param map a data frame or a matrix of dimension p and 2.
#'  The row names should identical to the row names of \code{haplo}.
#' @param RR a numeric value corresponding de the recombination rate.
#' @param h2 a numeric value [0,1] or vector of length y.
#' @param nb_gen an integer indicating the number of generation
#'   to simulate from the current population. Default = 15
#' @param nb_ind_gen an integer or a vector of the same length \code{nb_gen}
#'   indicating the number of individual per generation. Default = 2000
#' @param nb_ind_P0 an integer indicating the final number of individuals in the population
#' @param folds an integer for the number of folds. Default  = 10.
#' @param n_rep an integer for the number of replicates. Default  = 10.
#' @param method a character string indicating which method to use
#'   for the calibration of the genotype-phenotype relationship:
#'   "RF" (default), "SVR", or "NW".
#' @param nb_cores the number of cores to use for parallel execution. Default = 1.
#' @param verbose logical scalar: TRUE (default) or FALSE
#' @param ... further arguments to be passed to or from methods.
#'
#' @return A list :
#'  *
#'  *
#' @export
#'
#' @examples
simul_P0 <- function(haplo,
                     pheno,
                     map,
                     RR,
                     h2,
                     nb_gen = 15,
                     nb_ind_gen = 2000,
                     nb_ind_P0 = 500,
                     folds = 10,
                     n_rep = 10,
                     method = c('RF', 'SVR', 'NW'),
                     nb_cores = 1,
                     verbose = T,
                     ...) {
  ellipsis <- list(...)
  if (is.null(ellipsis$ntree)) {
    ntree <- 1000
  } else{
    ntree <- ellipsis$ntree
  }
  if (is.null(ellipsis$mtry)) {
    mtry <- floor(sqrt(nrow(haplo)))
  } else{
    mtry <- ellipsis$mtry
  }
  if (is.null(ellipsis$kernel_type)) {
    kernel_type <- 'aggregate'
  } else{
    kernel_type <- ellipsis$kernel_type
  }
  pheno <- as.matrix(pheno)
  method <- match.arg(method)
  geno <- haplo[, , 1] + haplo[, , 2]
  geno <- t(geno)
  if (verbose)
    cat("  Creating genotypes of the P0 population\n\n")
  P0 <- simul_geno(
    haplo = haplo,
    map = map,
    RR = RR,
    nb_gen = nb_gen,
    nb_ind_gen = nb_ind_gen,
    nb_ind_P0 = nb_ind_P0,
    nb_cores = nb_cores,
    verbose = verbose
  )
  geno_P0 <- P0$sim_geno

  nb_trait <- ncol(pheno)
  if (any(lengths(ellipsis) > 1)) {
    if (verbose)
      cat("\n  Optimization of parameters for all the traits\n\n")
    opt_param <- vector("list", nb_trait)
    for (tr in 1:nb_trait) {
      if (verbose) {
        cat("Trait ", tr, " over", nb_trait, '\n')
      }
      if (method == 'RF') {
        mod <- rf_opt(
          x = geno,
          y = pheno[, tr],
          ntree = ntree,
          mtry = mtry,
          folds = folds,
          n_rep = n_rep,
          nb_cores = nb_cores
        )
      } else if (method == 'SVR') {
        mod <- svr_opt(
          x = geno,
          y = pheno[, tr],
          C = ellipsis$C,
          sigma = ellipsis$sigma,
          folds = folds,
          n_rep = n_rep,
          nb_cores = nb_cores
        )
      } else if (method == 'NW') {
        mod = nw_opt(
          x = geno,
          y = as.vector(pheno[, tr]),
          kernel_type = kernel_type,
          h_grid = ellipsis$h,
          folds = folds,
          n_rep = n_rep,
          nb_cores = nb_cores
        )
      }
      opt_param[[tr]] <- mod$bestTune
    }
    if (verbose) {
      cat("\n  Optimized parameters are \n")
      print(do.call('rbind', opt_param))
    }
  } else {
    cat("\n  No parameters optimization \n")
    opt_param <- rep(list(ellipsis), nb_trait)
    print(do.call('rbind', opt_param))
  }

  if (verbose)
    cat("\n  Creating phenotypes of the P0 population using \n\n")
  gen_val_P0 <- vector("list", nb_trait)
  phen_val_P0 <- vector("list", nb_trait)
  sigma2_eps <- rep(0, nb_trait)
  for (tr in 1:nb_trait) {
    if (verbose) {
      cat("Trait ", tr, " over", nb_trait, '\n')
    }
    if (method == 'RF') {
      gen_val_P0[[tr]] <- rf(
        x = geno,
        y = pheno[, tr],
        x_pred = geno_P0,
        ntree = ntree,
        mtry = opt_param[[tr]]$mtry
      )
    } else if (method == 'SVR') {
      gen_val_P0[[tr]] <- svr(
        x = geno,
        y = pheno[, tr],
        x_pred = geno_P0,
        C = opt_param[[tr]]$C,
        sigma = opt_param[[tr]]$sigma
      )
    } else if (method == 'NW') {
      gen_val_P0[[tr]] <- nw(
        x = geno,
        y = as.vector(pheno[, tr]),
        x_pred = geno_P0,
        h = opt_param[[tr]]$h,
        kernel_type = kernel_type,
        opt_subset = NULL
      )
    }
    sigma2_eps[tr] <-
      var(gen_val_P0[[tr]]) * ((1 - h2[tr]) / h2[tr])
    phen_val_P0[[tr]] <- gen_val_P0[[tr]] + rnorm(nb_ind_P0,
                                                  mean = 0,
                                                  sd = sqrt(sigma2_eps[tr]))
  }
  return(
    list(
      "haplo" = P0$sim_haplo,
      "geno" = geno_P0,
      "gen_val" = gen_val_P0,
      "phen_val" = phen_val_P0,
      "map" = map,
      "opt_param" = opt_param
    )
  )
}
