# simul_P0 ----------------------------------------------------------------
#' Generate a reference population (P0) from real genotypic and phenotypic data.
#'
#' \code{simul_geno} compute...
#'
#' @param haplo an array of dimension p, n and 2 or a matrix of dimension p, 2n (NAs not  allowed for NW method).
#' @param map a data frame or a matrix of dimension p and 2.
#'  The row names should identical to the row names of \code{haplo}.
#' @param RR a numeric value corresponding de the recombination rate.
#' @param nb_gen an integer indicating the number of generation
#'   to simulate from the current population. Default = 15
#' @param nb_ind_gen an integer or a vector of the same length \code{nb_gen}
#'   indicating the number of individual per generation. Default = 2000
#' @param nb_ind_P0 an integer indicating the final number of individuals in the population
#' @param nb_cores the number of cores to use for parallel execution. Default = 1.
#' @param verbose logical scalar: TRUE (default) or FALSE
#'
#' @return A list with two object:
#'   * an array containing the haplotypes
#'   * a matrix with the genotypes
#' @export
#'
#' @examples
#'
#' # Takes time!
#' system.time(
#'   P0 <- simul_geno(
#'     haplo = rice1_x,
#'     map = rice1_map,
#'     RR = 244000,
#'     nb_ind_gen = 1000,
#'     nb_gen = 2,
#'     nb_ind_P0 = 200,
#'   )
#'  )
#'
#' dim(P0$sim_haplo)
#' dim(P0$sim_geno)
#'
#'
simul_geno <-
  function(haplo,
           map,
           RR,
           nb_gen = 15,
           nb_ind_gen = 2000,
           nb_ind_P0 = 500,
           nb_cores = 1,
           verbose = T) {
    n <- ncol(haplo)
    genome <- haplo
    rnames <- rownames(genome)
    if (length(nb_ind_gen) == 1) {
      nb_ind_gen <- rep(nb_ind_gen, nb_gen)
    }
    for (i in 1:nb_gen) {
      if (verbose) {
        cat("Generation", i, "over", nb_gen, '\n')
      }
      sim_haplo <-
        multi_cross(genome, map, RR, nb_ind_gen[i], nb_cores)
      genome <- sim_haplo
    }
    if (!is.null(nb_ind_P0)) {
      sim_haplo <-
        sim_haplo[, sample(ncol(sim_haplo), nb_ind_P0, replace = F), ]
    }
    rownames(sim_haplo) <- rnames
    sim_haplo <- aperm(sim_haplo, perm= c(2,1,3))
    sim_geno <- sim_haplo[, , 1] + sim_haplo[, , 2]
    return(list(sim_haplo = sim_haplo, sim_geno = sim_geno))
  }


# multi_cross -------------------------------------------------------------------

multi_cross <- function(genome, map, RR, nb_ind_gen, nb_cores) {
  ped <- matrix(0, nrow = nb_ind_gen, ncol = 2)
  rownames(ped) <- 1:nb_ind_gen
  for (m in 1:nb_ind_gen) {
    ped[m, ] <- sample(1:(dim(genome)[2]), size = 2, replace = FALSE)
  }
  if (nb_cores > 1) {
    doParallel::registerDoParallel(cores = nb_cores)
    '%dopar%' <- foreach::`%dopar%`
    off_haplo <-
      foreach::foreach(m = 1:nb_ind_gen) %dopar% {
        cross(ped[m, ], genome, map, RR)
      }
  } else {
    off_haplo <- replicate(n = nb_ind_gen,
                           expr = {
                             array(dim = c(dim(genome)[1], 1, dim(genome)[3]))
                           },
                           simplify = F)
    for (m in 1:nb_ind_gen) {
      off_haplo[[m]] <- cross(ped[m, ], genome, map, RR)
    }
  }
  off_haplo <- abind::abind(off_haplo, along = 2)
  return(off_haplo)
}


# cross -------------------------------------------------------------------

cross <- function(ped, genome, map, RR) {
  parent1 <- ped[1]
  parent2 <- ped[2]
  parent_genome <- genome[, c(parent1, parent2), ]

  p_vec <- c(0, diff(map$position / RR) / 100)
  p_vec[p_vec < 0] <- rbinom(sum(p_vec < 0), 1, 0.5)

  rec_gam1 <- cumsum(rbinom(length(p_vec), 1, p_vec))
  v_int <- rep(0:1, length.out = length(unique(rec_gam1)))
  rec_gam1 <- v_int[match(rec_gam1, unique(rec_gam1))]
  rec_gam2 <- cumsum(rbinom(length(p_vec), 1, p_vec))
  v_int <- rep(0:1, length.out = length(unique(rec_gam2)))
  rec_gam2 <- v_int[match(rec_gam2, unique(rec_gam2))]

  ref_strand1 <- rbinom(1, 1, 0.5) + 1
  alt_strand1 <- setdiff(c(1, 2), ref_strand1)
  ref_strand2 <- (rbinom(1, 1, 0.5) + 1)
  alt_strand2 <- setdiff(c(1, 2), ref_strand2)

  off_haplo1 <-
    ifelse(rec_gam1 == 0, parent_genome[, 1, ref_strand1],
           parent_genome[, 1, alt_strand1])
  off_haplo2 <-
    ifelse(rec_gam2 == 0, parent_genome[, 2, ref_strand2],
           parent_genome[, 2, alt_strand2])
  off_haplo <- array(dim = c(length(off_haplo1), 1, 2))
  off_haplo[, 1, 1] <- off_haplo1
  off_haplo[, 1, 2] <- off_haplo2
  return(off_haplo)
}
