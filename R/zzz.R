#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name methyldeconvstartup
NULL

.onLoad <- function(libname, pkgname){
  cli::cli_alert("checking methyldeconv environment and dependencies")
  
  # We ensure to have reticulate
  if (!dir.exists(reticulate::miniconda_path())) {
    message("Setting python version in miniconda to be 3.10")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.10)
    message("Setting up miniconda environment..")
    suppressMessages(reticulate::install_miniconda())
  }
  
  # We ensure to have the r-methyldeconv env
  if (!("r-methyldeconv" %in% reticulate::conda_list()$name)) {
    message("Create conda evironment 'r-methyldeconv' for MethAtlas...")
    reticulate::conda_create("r-methyldeconv", python_version = "3.10")
    message("Install all python dependencies...")
    reticulate::py_install(packages = c("numpy", "pandas", "scipy", "matplotlib") , envname = "r-methyldeconv",  method = "conda", pip = T)
  }
  
  paths <- reticulate::conda_list()
  path <- paths[paths$name == "r-methyldeconv", 2][[1]]
  
  Sys.setenv(RETICULATE_PYTHON = path)
  reticulate::use_miniconda(condaenv = "r-methyldeconv", required = TRUE)
  reticulate::py_config()
  reticulate::configure_environment(pkgname, force = TRUE)
}