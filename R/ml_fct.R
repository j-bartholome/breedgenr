# Optimisation of RF ------------------------------------------------------
#' Optimisation of randomForest parameters
#'
#' The \code{rf_opt} function performs a grid search to optimise
#'   the mtry parameter mtry for randomForest method.
#'
#'
#' @param x A data frame or a matrix of predictors (genotypes) of dimension n, p.
#' @param y A response vector (phenotype) of length n.
#' @param ntree Number of trees (integer). Default is 1,000
#' @param mtry an integer indicating the number of variables randomly sampled at each split.
#'   Default is the (rounded down) square root of the number variables.
#' @param folds an integer for the number of folds. Default = 10.
#' @param n_rep an integer for the number of replicates. Default = 10.
#' @param nb_cores the number of cores to use for parallel execution. Default = 1.
#'
#' @return The model
#'
#' @export
#'
#' @seealso
#'  \code{\link[caret]{train}}, \code{\link[ranger]{ranger}}, \code{\link[randomForest]{randomForest}}
#'
#' @examples
#'
#'
#'
rf_opt <- function(x,
                   y,
                   ntree = 1000,
                   mtry,
                   folds = 10,
                   n_rep = 10,
                   nb_cores = 1) {
  val_ctrl <- caret::trainControl(method = "repeatedcv",
                                  number = folds,
                                  repeats = n_rep)
  val_grid <-  expand.grid(mtry = mtry,
                           min.node.size = 5,
                           splitrule = 'variance')
  doParallel::registerDoParallel(cores = nb_cores)
  mod <-
    caret::train(
      x = x,
      y = y,
      method = "ranger",
      num.trees = ntree,
      trControl = val_ctrl,
      tuneGrid = val_grid
    )
  doParallel::stopImplicitCluster()
  return(mod)
}

# Random Forest -----------------------------------------------------------
#' Prediction using randomForest
#'
#' The \code{rf} function performs prediction using new genotypic data (x_pred)
#'
#'
#' @param x A data frame or a matrix of predictors (genotypes) of dimension n, p.
#' @param y A response vector (phenotype) of length n.
#' @param x_pred A data frame or a matrix (m, p) of the same predictors os \code{x} for a new set of observations.
#' @param ntree Number of trees (integer). Default is 1,000
#' @param mtry an integer indicating the number of variables randomly sampled at each split.
#'   Default is the (rounded down) square root of the number variables.
#'
#' @return A vector of predicted values of length m
#'
#' @export
#'
#' @seealso #'  \code{\link[ranger]{ranger}}, \code{\link[randomForest]{randomForest}}
#'
#' @examples
#'
#'
rf <- function(x,
               y,
               x_pred,
               ntree = 1000,
               mtry) {
  ctrl <- caret::trainControl(method = "none")
    mod <-
      caret::train(
        x = x,
        y = y,
        method = "ranger",
        num.trees = ntree,
        tuneGrid = data.frame(
          mtry = mtry,
          min.node.size = 5,
          splitrule = 'variance'
        ),
        trControl = ctrl
      )
  pred_val <- predict(mod, newdata = x_pred)
  return(pred_val)
}


# Optimisation of SVR -----------------------------------------------------
#' Optimisation of support vector regression with gaussian kernel
#'
#' The \code{svr_opt} function performs a grid search to optimise
#'   the C and sigma parameters for support vector regression.
#'   The \code{\link[caret]{train}} function is used.
#'
#' @param x A data frame or a matrix of predictors (genotypes) of dimension n, p.
#' @param y A response vector (phenotype) of length n.
#' @param C  a numeric value for the soft margin parameter
#' @param sigma  a numeric value indicating the inverse kernel width for the gaussian kernel
#' @param folds an integer for the number of folds. Default = 10.
#' @param n_rep an integer for the number of replicates. Default= 10.
#' @param nb_cores the number of cores to use for parallel execution. Default = 1.
#'
#' @return The model
#' @export
#'
#' @examples
svr_opt <- function(x,
                    y,
                    C,
                    sigma = NULL,
                    folds = 10,
                    n_rep = 10,
                    nb_cores = 1) {
  if (is.null(sigma)) {
    df <- data.frame(x, y)
    sigma <- unname(round(kernlab::sigest(y ~ ., data = df), 5))
  }
  val_ctrl <- caret::trainControl(method = "repeatedcv",
                                  number = folds,
                                  repeats = n_rep)
  val_grid <-  expand.grid(sigma = sigma,
                           C = C)
  doParallel::registerDoParallel(cores = nb_cores)
  mod <-
    caret::train(
      x = x,
      y = y,
      method = "svmRadial",
      trControl = val_ctrl,
      tuneGrid = val_grid
    )
  doParallel::stopImplicitCluster()
  return(mod)
}


# SVR ---------------------------------------------------------------------
#' Prediction using support vector regression with gaussian kernel
#'
#' The \code{svr} function performs prediction using new genotypic data (x_pred)
#'
#'
#' @param x A data frame or a matrix of predictors (genotypes) of dimension n, p.
#' @param y A response vector (phenotype) of length n.
#' @param x_pred A data frame or a matrix (m, p) of the same predictors os \code{x} for a new set of observations.
#' @param C  a numeric value for the soft margin parameter
#' @param sigma  a numeric value indicating the inverse kernel width for the gaussian kernel
#'
#' @return  A vector of predicted values of length m
#' @export
#'
#' @examples
svr <- function(x,
                y,
                x_pred,
                C,
                sigma) {
  val_ctrl <- caret::trainControl(method = "none")
  val_grid <-  expand.grid(sigma = sigma,
                           C = C)
  mod <-
    caret::train(
      x = x,
      y = y,
      method = "svmRadial",
      trControl = val_ctrl,
      tuneGrid = val_grid
    )
  pred_val <- predict(mod, newdata = x_pred)
  return(pred_val)
}
