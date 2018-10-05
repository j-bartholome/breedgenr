# Evaluation with package data  -------------------------------------------

# data prep
haplo <- rice1_x
x <- aperm(haplo, perm = c(2, 1, 3))
x <- x[, , 1] + x[, , 2]
y <- rice1_y[, 1]
pheno <- rice1_y
map <- rice1_map
RR <- 244000
h2 <- rep(0.7, 3)
nb_gen <- 5
nb_ind_gen <- 500
nb_ind_P0 <- 100
folds <- 5
n_rep <- 2
nb_cores <- 6
verbose <- T

dim(x)
dim(haplo)

# P0 genotype generation
system.time(
  P0 <- simul_geno(
    haplo = haplo,
    map = map,
    RR = RR,
    nb_ind_gen = nb_ind_gen,
    nb_gen = nb_gen,
    nb_ind_P0 = nb_ind_P0,
    nb_cores = nb_cores,
    verbose = T
  )
)

dim(P0$sim_haplo)
dim(P0$sim_geno)

# RF ----------------------------------------------------------------------
system.time(
  param_opt <- rf_opt(
    x = x,
    y = y,
    ntree = 100,
    mtry = c(floor(sqrt(nrow(rice1_x))^c(0.8, 1, 1.2))),
    folds = folds,
    n_rep = n_rep,
    nb_cores = nb_cores
  )
)
system.time(P0_val <- rf(
  x = x,
  y = y,
  x_pred = P0$sim_geno,
  ntree = 1000,
  mtry = param_opt$bestTune$mtry
))
system.time(
  P0_RF <- simul_P0(
    haplo = haplo,
    pheno = pheno,
    map = map,
    RR = RR,
    h2 = h2,
    nb_gen = nb_gen,
    nb_ind_gen = nb_ind_gen,
    nb_ind_P0 = nb_ind_P0,
    method = 'RF',
    ntree = 500,
    mtry = c(floor(sqrt(nrow(rice1_x))^c(0.8, 1, 1.2))),
    folds = folds,
    n_rep = n_rep,
    nb_cores = nb_cores,
    verbose = verbose
  )
)

plot(P0_RF$gen_val[[1]], P0_RF$phen_val[[1]])


# SVR ---------------------------------------------------------------------
system.time(
  param_opt <- svr_opt(
    x = x,
    y = y,
    C = 10 ^ (2),
    sigma = kernlab::sigest(y ~ x),
    folds = folds,
    n_rep = n_rep,
    nb_cores = nb_cores
  )
)
system.time(
  P0_val <- svr(
    x = x,
    y = y,
    x_pred = P0$sim_geno,
    C = param_opt$bestTune$C,
    sigma = param_opt$bestTune$sigma
  )
)
system.time(
  P0_SVR <- simul_P0(
    haplo = haplo,
    pheno = pheno,
    map = map,
    RR = RR,
    h2 = h2,
    nb_gen = nb_gen,
    nb_ind_gen = nb_ind_gen,
    nb_ind_P0 = nb_ind_P0,
    method = 'SVR',
    C = 10 ^ (1:3),
    sigma = kernlab::sigest(y ~ x),
    folds = folds,
    n_rep = n_rep,
    nb_cores = nb_cores,
    verbose = verbose
  )
)

plot(P0_SVR$gen_val[[1]], P0_SVR$phen_val[[1]])

# NW ----------------------------------------------------------------------
system.time(
  h_opt <- nw_opt(
    x = x,
    y = y,
    kernel_type = 'aggregate',
    h_grid = kernlab::sigest(y ~ x, frac = 1, scaled = F),
    n_rep = n_rep,
    folds = folds,
    nb_cores = nb_cores
  )
)

system.time(pred_NW <- nw(
  x = x,
  y = y,
  h = 0.1,
  x_pred = P0$sim_geno
))

system.time(
  P0_NW <- simul_P0(
    haplo = haplo,
    pheno = pheno,
    map = map,
    RR = RR,
    h2 = h2,
    nb_gen = nb_gen,
    nb_ind_gen = nb_ind_gen,
    nb_ind_P0 = nb_ind_P0,
    method = 'NW',
    kernel_type = kernel_type,
    h = kernlab::sigest(x),
    n_rep = n_rep,
    folds = folds,
    nb_cores = nb_cores,
    verbose = verbose
  )
)

plot(P0_NW$gen_val[[1]], P0_NW$phen_val[[1]])



# Parameter evaluation ----------------------------------------------------

pop <- create_population(haplo = haplo,
                         pheno = NULL)

pop <- create_population(haplo = haplo,
                         pheno = pheno)

pop <- create_population(haplo = haplo,
                         pheno = pheno,
                         map = map)

pop <- create_population(haplo = haplo,
                         pheno = pheno,
                         map = map,
                         h2 = 0.3)

pop <- create_population(
  haplo = haplo,
  pheno = pheno,
  map = map,
  h2 = rep(0.3, 3),
  RR = RR
)

pop <- create_population(
  haplo = haplo,
  pheno = pheno,
  map = map,
  h2 = rep(0.3, 3),
  RR = RR,
  mtry = 10
)

pop <- create_population(
  haplo = haplo,
  pheno = pheno,
  map = map,
  h2 = rep(0.3, 3),
  RR = RR,
  mtry = 10
)

pop <- create_population(
  haplo = haplo,
  pheno = pheno,
  map = map,
  RR = RR,
  h2 =rep(0.3, 3),
  nb_gen = 2,
  nb_ind_gen = 200,
  nb_ind_P0 = 100,
  mtry = 40
)


pop <- create_population(
  haplo = haplo,
  pheno = pheno,
  map = map,
  RR = RR,
  h2 = 0.4,
  simul = TRUE,
  nb_gen = 2,
  nb_ind_gen = 200,
  nb_ind_P0 = 100,
  folds = 5,
  n_rep = 2,
  nb_cores = 4,
  verbose = T,
  mtry = 40
)


pop <- create_population(
  haplo = haplo,
  pheno = pheno,
  map = map,
  RR = RR,
  h2 = 0.4,
  simul = TRUE,
  method = 'TT',
  nb_gen = 15,
  nb_ind_gen = 2000,
  nb_ind_P0 = 500,
  folds = 5,
  n_rep = 5,
  nb_cores = 1,
  verbose = T
)
