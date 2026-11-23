# load_engine.R -- compile (first call only) and load cpp/sbart_engine.cpp.
# The compiled library is kept in cpp/build (or in $SBART_ENGINE_CACHE), so jobs
# only load it. Run setup.sh once before submitting jobs, so that the jobs do not
# all compile at the same time.
load_engine <- function(code_dir, rebuild = FALSE, verbose = FALSE) {
  cpp   <- file.path(code_dir, "cpp", "sbart_engine.cpp")
  cache <- Sys.getenv("SBART_ENGINE_CACHE", file.path(code_dir, "cpp", "build"))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  # std::thread needs -pthread on older Linux systems
  old <- Sys.getenv("PKG_LIBS", unset = NA)
  Sys.setenv(PKG_LIBS = trimws(paste(if (is.na(old)) "" else old, "-pthread")))
  on.exit(if (is.na(old)) Sys.unsetenv("PKG_LIBS") else Sys.setenv(PKG_LIBS = old))
  # sourceCpp rewrites the cache index (cache.rds) at every load; a job that reads it
  # while another job writes it gets an error, so try again after a short pause
  for (k in 1:10) {
    ok <- tryCatch({ Rcpp::sourceCpp(cpp, cacheDir = cache, env = globalenv(), rebuild = rebuild,
                                     showOutput = verbose, verbose = verbose); TRUE },
                   error = function(e) { if (k == 10) stop(e); FALSE })
    if (ok) break
    Sys.sleep(runif(1, 1, 5))
  }
  invisible(TRUE)
}
