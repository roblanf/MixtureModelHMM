#' Generates report
#'
#' @param hmm_result predict_tree/predict_tree_mixed returned list
#' @param output_filename name of file to be saved
#'
#' @return saves prediction report in .txt format
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' save_report(hmm_result,"report")
save_report<-function(hmm_result,output_filename){
  tained_model<-hmm_result[[3]]
  sink(paste(output_filename,'.txt',sep = ""))
  cat("=============================\n")
  cat("Probability of transition of one class to another\n")
  cat("Each class is represented by C\n")
  print(tained_model$A)
  cat("=============================\n")
  cat("Probability of emission from each class\n")
  cat("Each emission is represented by E\n")
  cat("Emissions depend on the model selected\n")
  cat("Model 1 has no additional emissions\n")
  cat("Model 2 has one additional emissions for constant sites\n")
  cat("Model 3 has two additional emissions for constant sites and non-informative sites\n")
  cat("Model 4 has three additional emissions for constant sites, non-informative sites and same parsimony sites\n")
  cat("Model 5 has additional emissions depending which model(from 1-4) it converged on\n")
  cat("Model 4 does not apply for mixture models as all classes have same parsimony\n")
  print(tained_model$E)
  cat("=============================\n")
  classification<- hmm_result$classification
  ###############################
  v.pred<- hmm_result[[1]];p.pred<-hmm_result[[2]]
  cat("Predicted sites for each class\n")
  cat("Viterbi dynamic Program\n")
  print(table(v.pred))
  cat("Posterior decoding Program\n")
  print(table(p.pred))
  cat("=============================\n")
  cat("Proportion of predicted sites for each class\n")
  cat("Viterbi dynamic Program\n")
  print(prop.table(table(v.pred)))
  cat("Posterior decoding Program\n")
  print(prop.table(table(p.pred)))
  sink()
}
