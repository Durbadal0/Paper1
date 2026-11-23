# compare_outputs.R DIR_R DIR_CPP -- compares the replicate files of the R version
# (DIR_R) and of the C++ version (DIR_CPP) field by field. The run time is ignored.
args <- commandArgs(TRUE)
fr <- sort(list.files(args[1], "^rep[0-9]+\\.rds$"))
if (!length(fr)) stop("no replicate files in ", args[1])
worst <- 0; n_ident <- 0
for (f in fr) {
  a <- readRDS(file.path(args[1], f)); b <- readRDS(file.path(args[2], f))
  a$minutes <- b$minutes <- NULL
  d <- sapply(names(a), function(n) {
    if (identical(a[[n]], b[[n]])) return(0)
    x <- suppressWarnings(as.numeric(unlist(a[[n]]))); y <- suppressWarnings(as.numeric(unlist(b[[n]])))
    if (length(x) != length(y) || !identical(is.na(x), is.na(y))) return(Inf)
    max(abs(x - y), na.rm = TRUE) })
  same <- all(d == 0); n_ident <- n_ident + same
  worst <- max(worst, d)
  cat(sprintf("%s  scenario %s: %s\n", f, a$scenario,
              if (same) "all fields identical" else
              paste(sprintf("%s differs by up to %.2g", names(d)[d > 0], d[d > 0]), collapse = "; ")))
  cat(sprintf("        AMSE  R: %s\n        AMSE C++: %s\n",
              paste(sprintf("%.10f", a$amse), collapse = " "), paste(sprintf("%.10f", b$amse), collapse = " ")))
}
cat(sprintf("\n%d of %d replicates identical in every field; largest difference in any field: %.2g\n",
            n_ident, length(fr), worst))
if (!is.finite(worst) || worst > 1e-10) quit(status = 1)
