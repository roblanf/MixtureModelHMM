#' Generates report
#'
#' @param hmm_result run_hmm() returned object
#' @param output_filename name of file to be saved
#'
#' @return saves prediction report in HTML file
#' @importFrom rmarkdown draft render
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' save_report(hmm_result,"hmm_report.Rmd")
save_report<-function(hmm_result, output_filename){
  draft(output_filename,template = "HMM_report", package="MixtureModelHMM")
  render(output_filename,params = list(hmm_result=hmm_result))
}
