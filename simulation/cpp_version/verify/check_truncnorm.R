# check_truncnorm.R -- the engine's truncated-normal sampler gives the same draws as
# truncnorm::rtruncnorm for the two cases the sampler uses: (0, Inf) and (-Inf, 0).
#   Rscript verify/check_truncnorm.R        (run from the folder that holds code/)
suppressMessages(library(truncnorm))
source("code/cpp/load_engine.R"); load_engine("code")
set.seed(5)
mu <- c(runif(2e5, -6, 6), rnorm(1e5, 0, 0.3), seq(-0.46, -0.44, length.out = 1000),
        seq(0.44, 0.46, length.out = 1000))
set.seed(9); a1 <- rtruncnorm(length(mu), a = 0, b = Inf, mean = mu, sd = 1)
             a2 <- rtruncnorm(length(mu), a = -Inf, b = 0, mean = mu, sd = 1); ua <- runif(1)
set.seed(9); b1 <- engine_rtruncnorm(mu, 1L); b2 <- engine_rtruncnorm(mu, 0L); ub <- runif(1)
ok <- identical(a1, b1) && identical(a2, b2) && identical(ua, ub)
cat(sprintf("truncated normal, %d draws on each side: %s\n", length(mu),
            if (ok) "IDENTICAL to truncnorm::rtruncnorm" else "DIFFERENT"))
if (!ok) quit(status = 1)
