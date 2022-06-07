#' Saves the predicted output to gz file
#'
#' @param hmm_result predict_tree/predict_tree_mixed returned list
#' @param file_name name of file to be saved
#' @param algo "viterbi" or "posterior"
#'
#' @return saves prediction output in .gz file
#' @export
#'
#' @examples
#' save_file(pred,"viterbi","output_v")
#' save_file(pred,"posterior","output_p")
#' v.pred<-system("gzcat output_v.gz",intern=TRUE)
save_file<-function(hmm_result,algo,file_name){

  if(missing(algo)) {
    algo="viterbi"
  }
  if(algo=="viterbi"){
    vit.pred=hmm_result[[1]]
    df1=data.frame(vit.pred)
  }
  else if(algo=="posterior"){
    post.pred=hmm_result[[2]]
    df1=data.frame(post.pred)
  }
  else{
    print("Invalid algorithm choice")
  }
  gz1 <- gzfile(paste(file_name,".gz",sep=''), "w")
  write.table(df1, gz1,row.names = F, col.names = F,quote=F)
  close(gz1)

}


#' Generates report
#'
#' @param hmm_result predict_tree/predict_tree_mixed returned list
#' @param output_filename name of file to be saved
#'
#' @return saves prediction report in .txt
#' @export
#'
#' @examples
#' save_report(hmm_result,"report")
save_report<-function(hmm_result,output_filename){
  tained_model<-hmm_result[[5]]
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
