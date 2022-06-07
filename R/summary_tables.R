#' Create a transition table
#'
#' Where each row has:
#' *site
#' *class_to
#' *class_from
#' And we ignore the first site, because it's not a transition.
#'
#' @param hmm_result run_hmm() returned object
#'
#' @return dataframe of transition table
#' @export
#'
#' @examples
#' transition_table(hmm_result)
transition_table<-function(hmm_result){
  site=head(cumsum(rle(hmm_result$classification)$lengths)+1, -1)
  class_from=head(rle(hmm_result$classification)$values, -1)
  class_to=tail(rle(hmm_result$classification)$values, -1)
  data.frame(site,class_from,class_to)
}

#' Generate partition file
#'
#' @param hmm_result run_hmm() returned object
#' @param output_filename name of file to be saved
#'
#' @return saves partition file in .nex format
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' partition_file(hmm_result,"hmm_partitions")
partition_file<-function(hmm_result,output_filename){
  classification<-hmm_result$classification
  site=cumsum(rle(hmm_result$classification)$lengths)
  site=c(1,site)
  from_site=head(site,-1)
  to_site=tail(site,-1)
  classes=as.numeric(gsub("C","",rle(hmm_result$classification)$values))
  partition_table=data.frame(from_site,to_site,classes)
  partition_table$from_to_site=paste(partition_table$from_site,partition_table$to_site,sep='-')
  partition_tibble<-partition_table %>% group_by(classes) %>% summarise(lst= paste0(from_to_site, collapse = " "))
  sink(paste(output_filename,'.nex',sep = ""))
  cat("#nexus\n")
  cat("begin sets;\n")
  for(i in 1:length(unique(classification))){
    cat(paste("\tcharset class_",i," = ",partition_tibble[i,]$lst,";\n",sep = ""))
  }
  cat("end;\n")
  sink()
}


#' Summary
#'
#' @param hmm_result run_hmm() returned object
#' @param ... additional arguments to be passed
#'
#' @return summary
#' @export
#'
#' @examples
#' summary(hmm_result)
summary.MixtureModelHMM<-function(hmm_result, ...){

  print(paste("Number of sites: ",nrow(hmm_result$data)))
  print(paste("Number of classes: ",ncol(hmm_result$data)))
  print("Sites in each class: ")
  print(table(hmm_result$classification))
  print(paste("Number of transitions: ",length(rle(hmm_result$classification)$value)-1))
  print(paste("Algorithm used to determine the sequence: ",hmm_result$algorithm))
  print(paste("Input files: ",hmm_result$site_input_file,hmm_result$aln_input_file))

}

