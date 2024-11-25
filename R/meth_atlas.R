#' run meth_atlas
#'
#' @param methyl_set A minfi MethylSet
#' @param reference_atlas Path to a csv file that saves a reference matrix with CpGs as rows and cell types as columns.
#'                        The default reference file is stored in 'inst/reference_atlas.csv'. 
#' @param temp_dir Path to directory where the beta matrix will be saved as a csv file.
#' @param out_dir Path to output directory. Output will be a csv file and a png representing the cell type fractions.
#' 
#' @export
#'
run_meth_atlas <- function(beta_matrix, reference_atlas = system.file("reference_atlas.csv", package = "methylDeconv"), temp_dir = NULL, out_dir = NULL){
  # check if python is installed, else install
  init_python()
  
  # set up temporary nd output directories
  tmp_dir <- temp_dir
  if (is.null(temp_dir)) {
    tmp_dir <- tempdir()
    dir.create(tmp_dir, showWarnings = FALSE)
  }
  
  if (is.null(out_dir)) {
    out_dir <- tmp_dir
  } 
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE)
  }
  
  # create a beta matrix from the methyl_set and save to temporary folder
  beta_matrix <- check_input_beta(beta_matrix)
  beta_path = paste0(tmp_dir, "/beta.csv")
  write.csv(beta_matrix, beta_path)
  
  # run meth_atlas
  system(paste("python", system.file("deconvolve.py", package = "methylDeconv")," -a", reference_atlas, beta_path, "--out", out_dir))
  
  # read the results to provide as data frame
  as.matrix(t(read.csv(paste0(out_dir, "/beta_deconv_output.csv"),
                      row.names = 1, check.names = FALSE
  )))
}
  