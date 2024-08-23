#' run meth_atlas
#'
#' @param meth methylated data matrix
#' @param unmeth unmethylated data matrix
#' @param reference_atlas Path to a csv file that saves a reference matrix with CpGs as rows and cell types as columns.
#'                        The default reference file is stored in 'inst/reference_atlas.csv'. 
#' @param temp_dir Path to directory where the beta matrix will be saved as a csv file.
#' @param out_dir Path to output directory. Output will be a csv file and a png representing the cell type fractions.
#' 
#' @return
#' @export
#'
#' @examples
run_meth_atlas <- function(meth, unmeth, reference_atlas = system.file("reference_atlas.csv", package = "methylDeconv"), temp_dir = NULL, out_dir = NULL){
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
  
  mSet <- minfi::MethylSet(meth, unmeth)
  
  # create a beta matrix from the mSet and save to temporary folder
  beta <- minfi::getBeta(mSet)
  beta_path = paste0(tmp_dir, "/beta.csv")
  write.csv(beta, beta_path)
  
  # run meth_atlas
  system(paste("python", system.file("deconvolve.py", package = "methylDeconv")," -a", reference_atlas, beta_path, "--out", out_dir))
  
  # read the results to provide as data frame
  t(utils::read.table(paste0(out_dir, "/beta_deconv_output.csv"),
                      sep = ",", header = TRUE,
                      row.names = 1, check.names = FALSE
  ))
}
  