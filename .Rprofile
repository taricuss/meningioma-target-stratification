local_lib <- normalizePath(file.path(getwd(), "r_lib"), winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}
rm(local_lib)
