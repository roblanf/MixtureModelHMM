#' Generates report
#'
#' Rmd file is generated which is rendered into HTML file
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

  #check for correct file extension
  if(unlist(strsplit(output_filename,"[.]"))[2] == "Rmd"){

    # create a Rmd file from template inside the package
    draft(output_filename,template = "hmm_report", package="MixtureModelHMM")

    # convert to HTML file
    render(output_filename,params = list(hmm_result=hmm_result))
  }
  else{
    print("Please use file extension .Rmd")
  }
}
