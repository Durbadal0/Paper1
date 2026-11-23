# install_packages.R -- install the R packages the simulation needs, if missing.
# Packages go to the first library on .libPaths() (R_LIBS_USER_DIR in config.sh).
need <- c("SoftBart", "MASS", "truncnorm", "survival", "ranger")
lib  <- .libPaths()[1]
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  cat("installing:", miss, "into", lib, "\n")
  install.packages(miss, lib = lib, repos = "https://cloud.r-project.org")
}
for (p in need) cat(sprintf("%-10s %s\n", p, as.character(packageVersion(p))))
cat("R", R.version$major, ".", R.version$minor, "\n", sep = "")
cat("The paper's local check used R 4.6.1 with SoftBart 1.0.3, truncnorm 1.0.9, ranger 0.18.0,\n",
    "survival 3.8.6 and MASS 7.3.65. Other versions run the same code but may not give\n",
    "the same random draws.\n", sep = "")
