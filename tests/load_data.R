
# # Load  -------------------------------------------------------------------
# pheno <-
#   read.table(file = "./data/Phenotypes.txt", sep = "\t", header = T)
# rownames(pheno) <- pheno$GID.
# geno <-
#   read.table(file = "./data/Genotypes.txt", sep = "\t", header = T)
# geno <- as.matrix(geno[, -1])
# haplo <-
#   read.table(file = "./data/Phased_genotypes.txt", sep = " ", header = T)
# haplo <- haplo[, -c(1:2)]
# haplo <- as.matrix(haplo)
# map <-
#   read.table(file = "./data/Physical_map.txt", sep = " ", header = T)
# colnames(map) <- c("ID", "chromosome", "position")
# map_centro <-
#   read.table(file = "./data/Physical_map_centromeres.txt", sep = " ", header = T)
#


# Evaluation with BGLR Data -----------------------------------------------
# data prep
data(wheat, package = 'BGLR')
x <- wheat.X
rownames(x) <- 1:nrow(x)
X <- t(x)
y <- wheat.Y[, 1]
haplo <- array(c(X, X), dim = c(dim(X)[1], dim(X)[2], 2),
               dimnames = list(rownames(X), colnames(X), NULL))
chr <- runif(nrow(haplo), min = 1, max = 12) %/% 1
mb <- runif(nrow(haplo), min = 0, max = 40000) %/% 1
map <- data.frame(chromosome = chr, position = mb)
map <- map[order(map$chromosome), ]
map$position <-
  do.call(c, tapply(map$position, map$chromosome, FUN = cumsum))
h2 <- rep(0.5, 4)
nb_ind_gen <- 200
nb_gen <- 10
nb_cores <- 8
RR <- 244
nb_ind_P0 <- 100
verbose <- T

dim(x)
dim(X)
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

# Kernel NW
system.time(
  h_opt <- nw_opt(
    x = x,
    y = y,
    kernel_type = 3,
    h_grid = seq(0.01, 1, by = 0.5),
    n_rep = 4,
    folds = 3,
    nb_cores = 8
  )
)

system.time(
  pred_NW <- nw(
    x = x,
    y = y,
    h = 0.1,
    x_pred = P0$sim_geno
  )
)

system.time(
  P0_NW <- simul_P0(
    haplo = haplo,
    pheno = y,
    map = map,
    RR = 100,
    h2 = 0.4,
    nb_gen = 10,
    nb_ind_gen = 500,
    nb_ind_P0 = 100,
    method = 'NW',
    kernel_type = 3,
    h = seq(0.01, 1, by = 0.5),
    n_rep = 4,
    folds = 3,
    nb_cores = 8,
    verbose = T
  )
)

plot(P0_NW$gen_val[[1]], P0_NW$phen_val[[1]])

# evaluation ML - RF
system.time(param_opt <- rf_opt(
  x = x,
  y = y,
  ntree = 100,
  mtry = seq(10, 100, 50),
  folds = 5,
  n_rep = 10,
  nb_cores = 8
))
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
    pheno = y,
    map = map,
    RR = 100,
    h2 = 0.4,
    nb_gen = 10,
    nb_ind_gen = 500,
    nb_ind_P0 = 100,
    method = 'RF',
    ntree = 500,
    mtry = seq(10, 100, 30),
    folds = 5,
    n_rep = 5,
    nb_cores = 8,
    verbose = T
  )
)

plot(P0_RF$gen_val[[1]], P0_RF$phen_val[[1]])

# evaluation ML - SVR
system.time(param_opt <- svr_opt(
  x = x,
  y = y,
  C = 10^(2),
  sigma = kernlab::sigest(y~x),
  folds = 5,
  n_rep = 10,
  nb_cores = 8
))
system.time(P0_val <- svr(
  x = x,
  y = y,
  x_pred = P0$sim_geno,
  C = param_opt$bestTune$C,
  sigma = param_opt$bestTune$sigma
))
system.time(
  P0_SVR <- simul_P0(
    haplo = haplo,
    pheno = y,
    map = map,
    RR = 100,
    h2 = 0.4,
    nb_gen = 10,
    nb_ind_gen = 500,
    nb_ind_P0 = 100,
    method = 'SVR',
    C = 10^(1:3),
    sigma = kernlab::sigest(y~x),
    folds = 5,
    n_rep = 5,
    nb_cores = 8,
    verbose = T
  )
)

plot(P0_SVR$gen_val[[1]], P0_SVR$phen_val[[1]])

test <- function(...) {
  ellipsis <- list(...)
  print(ellipsis$ntree)
  print(ellipsis)
  print(ellipsis$foo)
  any(lengths(ellipsis) > 0)
}

test(foo = 'a', ntree = 1000)

