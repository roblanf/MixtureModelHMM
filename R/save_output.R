#' Saves the predicted output to gz file
#'
#' @param pred predict_tree/predict_tree_mixed returned object
#' @param file_name name of file to be saved
#' @param algo "viterbi" or "posterior"
#'
#' @return saves prediction results
#' @export
#'
#' @examples
#' save_file(pred,"viterbi","output")
#' save_file(pred,"posterior","output")
#' v.pred<-system("gzcat viterbi.gz",intern=TRUE)
save_file<-function(pred,algo,file_name){

  if(missing(algo)) {
    algo="viterbi"
  }
  if(algo=="viterbi"){
    vit.pred=pred[[1]]
    df1=data.frame(vit.pred)
  }
  else if(algo=="posterior"){
    post.pred=pred[[2]]
    df1=data.frame(post.pred)
  }
  else{
    print("Invalid algorithm choice")
  }
  gz1 <- gzfile(paste(file_name,".gz",sep=''), "w")
  write.table(df1, gz1,row.names = F, col.names = F,quote=F)
  close(gz1)

}
