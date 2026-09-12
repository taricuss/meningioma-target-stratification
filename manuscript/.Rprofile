project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
local_lib <- file.path(project_root, "r_lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}
rm(project_root, local_lib)
