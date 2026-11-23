# check_forest_identity.R -- the engine's copy of the SoftBart forest code makes the
# same draws as the installed SoftBart package: same seed, 300 Gibbs sweeps, then the
# predictions, the split counts and the next random number are compared.
#   Rscript verify/check_forest_identity.R      (run from the folder that holds code/)
suppressMessages({library(SoftBart); library(Rcpp)})
Sys.setenv(PKG_LIBS = trimws(paste(Sys.getenv("PKG_LIBS"), "-pthread")))
sourceCpp("verify/forest_module.cpp")
set.seed(3); n <- 400; p <- 17
X <- matrix(runif(n * p), n, p); Y <- sin(3 * X[, 1]) + X[, 2] * X[, 3] + rnorm(n)
hyp <- Hypers(X, Y, num_tree = 200, sigma_hat = 1, k = 1/3)   # as in sim_worker.R
op  <- Opts(num_burn = 0, num_save = 1, update_sigma = FALSE)  # as in sim_worker.R
op2 <- op; op2$update_num_tree <- FALSE
Xt  <- matrix(runif(500 * p), 500, p)
run <- function(f) { set.seed(11)
  for (i in 1:300) { Z <- Y + rnorm(n, 0, 0.1); f$do_gibbs(X, Z, X, 1) }
  list(pred = f$do_predict(Xt), counts = as.numeric(f$get_counts()), u = runif(1)) }
a <- run(MakeForest(hyp, op, warn = FALSE))
b <- run(new(EngineForest, hyp, op2))
ok <- identical(a, b)
cat(sprintf("forest, 200 trees, 300 sweeps: %s (max |prediction difference| %g)\n",
            if (ok) "IDENTICAL to SoftBart" else "DIFFERENT", max(abs(a$pred - b$pred))))
if (!ok) quit(status = 1)
