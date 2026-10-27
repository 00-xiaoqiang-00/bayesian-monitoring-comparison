project_source_path <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
  starts <- c(if (nzchar(script) && file.exists(script)) dirname(normalizePath(script, mustWork = TRUE)) else character(), getwd())
  repository <- NULL
  for (start in unique(starts)) {
    root <- normalizePath(start, mustWork = TRUE)
    repeat {
      if (dir.exists(file.path(root, "code", "common")) && dir.exists(file.path(root, "data", "configuration"))) {
        repository <- root
        break
      }
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (!is.null(repository)) break
  }
  if (is.null(repository)) stop("Cannot locate the analysis repository.")
  function(...) {
    path <- file.path(repository, "code", "common", ...)
    if (!file.exists(path)) stop("Missing shared source: ", path)
    normalizePath(path, mustWork = TRUE)
  }
})
