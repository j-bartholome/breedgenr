# nw_optimization ---------------------------------------------------------
#' Optimisation of Nadaraya Watson regression parameters
#'
#' The \code{nw_opt} function performs a grid search to optimise
#'   the bandwith parameter for kernel estimation
#'
#'
#' @param x A data frame or a matrix of predictors (genotypes) of dimension n, p.
#' @param y A response vector (phenotype) of length n.
#' @param kernel_type a character string indicating which method to use
#'   for the kernel estimation: 'aggregate' (default), 'rho', 'gaussian'.
#' @param h_grid
#' @param folds an integer for the number of folds. Default = 10.
#' @param n_rep an integer for the number of replicates. Default = 10.
#' @param nb_cores the number of cores to use for parallel execution. Default = 1.
#'
#' @return A list
#' @export
#'
#' @examples
#' # May take some time
#' haplo <- rice1_x
#' x <- aperm(haplo, perm = c(2, 1, 3))
#' x <- x[, , 1] + x[, , 2]
#' y <- rice1_y[, 1]
#' system.time(
#'  h_opt <- nw_opt(
#'    x = x,
#'    y = y,
#'    kernel_type = 'aggregate',
#'    h_grid = kernlab::sigest(y ~ x, frac = 1),
#'    n_rep = 2,
#'    folds = 5,
#'    nb_cores = 2
#'  )
#'  )
#'
nw_opt <- function(x,
                   y,
                   kernel_type = c('rho', 'gaussian', 'aggregate'),
                   h_grid,
                   folds = 10,
                   n_rep = 10,
                   nb_cores = 1)
{
  kernel_type <- match.arg(kernel_type)
  n_ind <- length(y)
  v_mean_acc <- rep(0, length(h_grid))
  l <- 1
  n_rep <- n_rep*folds
  for (h  in h_grid)
  {
    cat("Evaluating the ",
        l,
        "occurence of h over",
        length(h_grid),
        '\n')
    cat('h value :', h, '\n')
    acc <- rep(0, n_rep)
    if (nb_cores > 1) {
      doParallel::registerDoParallel(cores = nb_cores)
      '%dopar%' <- foreach::`%dopar%`
      acc <-
        foreach::foreach(i = 1:n_rep, .combine = c) %dopar% {
          set.seed(i)
          validation_set <- sort(sample(n_ind,
                                        size = n_ind %/% folds,
                                        replace = FALSE))
          training_set <- c(1:n_ind)[-validation_set]
          y_val <- y[validation_set]
          x_val <- x[validation_set, ]
          y_pred <- kernel_nw(x[training_set, ],
                              y[training_set],
                              x_val,
                              h,
                              kernel_type)
          cor(y_val, y_pred)
        }
    } else {
      for (i in 1:n_rep)
      {
        set.seed(i)
        validation_set <- sort(sample(n_ind,
                                      size = n_ind %/% 5,
                                      replace = FALSE))
        training_set <- c(1:n_ind)[-validation_set]
        y_val <- y[validation_set]
        x_val <- x[validation_set, ]
        y_pred <- kernel_nw(x[training_set, ],
                            y[training_set],
                            x_val,
                            h,
                            kernel_type)
        acc[i] <- cor(y_val, y_pred)
      }
    }
    v_mean_acc[l] <- mean(acc, na.rm = T)
    if (l == 1) {
      max_temp <- v_mean_acc[l]
    }
    val_temp <- v_mean_acc[l]
    if (val_temp >= max_temp)
    {
      max_temp <- v_mean_acc[l]
      best_acc <- acc
    }
    l <- l + 1
  }
  return(list(
    bestTune = list(
      "h" = h_grid[which.max(v_mean_acc)],
      "opt_h_acc" = max(v_mean_acc),
      "opt_h_acc_range" = c(min(best_acc),
                            max(best_acc))
    )
  ))
}


# Kernel regression nw ----------------------------------------------------
#' Title
#'
#' @param x A data frame or a matrix of predictors (genotypes) of dimension n, p.
#' @param y A response vector (phenotype) of length n.
#' @param x_pred A data frame or a matrix (m, p) of the same predictors os \code{x} for a new set of observations.
#' @param h bandwith
#' @param kernel_type a character string indicating which method to use
#'   for the kernel estimation: 'aggregate' (default), 'rho', 'gaussian'.
#' @param opt_subset not used
#'
#' @return  A vector of predicted values of length m
#' @export
#'
#' @examples
#'
#'
#'
nw <-
  function(x,
           y,
           x_pred,
           h,
           kernel_type = c('rho', 'gaussian', 'aggregate'),
           opt_subset = NULL) {
    kernel_type <- match.arg(kernel_type)
    if (kernel_type == 'rho') {
      Kh_xi_xj <- rho_kernel
    } else if (kernel_type == 'gaussian') {
      Kh_xi_xj <- gaussian_kernel
    } else {
      Kh_xi_xj <- agg_kernel
    }
    if (is.vector(x_pred)) {
      x_pred <- as.matrix(x_pred, ncol = 1)
    }
    if (is.vector(x)) {
      x <- as.matrix(x, ncol = 1)
    }
    if (!is.null(opt_subset)) {
      x <- x[opt_subset,]
      y <- y[opt_subset]
    }
    mat_f <- mat_w <- matrix(nrow = nrow(x), ncol = nrow(x_pred))
    for (l in 1:nrow(x)) {
      mat_f[l, ] <- apply(x_pred, 1, Kh_xi_xj, xi = x[l,], h = h)
      mat_w[l, ] <- mat_f[l, ] * y[l]
    }
    norm_factor <- colSums(mat_f)
    weighted_pre <- colSums(mat_w)
    pred_val <- weighted_pre / norm_factor
    return(pred_val * 100)
  }

# rho_kernel --------------------------------------------------------------

rho_kernel <- function(xi, xj, h) {
  rho_kernel <- h * length(xi) * cor(xi, xj)
  return(rho_kernel)
}

# gaussian_kernel ---------------------------------------------------------

gaussian_kernel <- function(xi, xj, h) {
  gaussian_kernel <- exp(-h * (sum((xi - xj) ^ 2)))
  return(gaussian_kernel)
}

# agg_kernel --------------------------------------------------------------

agg_kernel <- function(xi, xj, h) {
  rho_kernel <- h * cor(xi, xj)
  gaussian_kernel <- exp(-h * (sum((xi - xj) ^ 2)))
  exph_kernel <- exp(-h * (sum(abs(xi - xj))))
  agg_kernel <- (rho_kernel + gaussian_kernel + exph_kernel)
  return(agg_kernel)
}
