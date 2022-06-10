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

  #RLE is Run length encoding to count each repeated element in the vector
  run_length=rle(hmm_result$classification)

  #Cummulative sum of runnning lengths to get the site number where each class is beginning from
  # we don't need the last elements as its the last site
  site=head(cumsum(run_length$lengths)+1, -1)

  #Get list of class transitions
  #The first class to second last class in list of run_length values are `from` classes
  class_from=head(run_length$values, -1)
  #The second class to last class in list of run_length values are `to` classes
  class_to=tail(run_length$values, -1)

  #return a dataframe
  data.frame(site,class_from,class_to)
}

#' Generate partition file
#'
#' @param hmm_result run_hmm() returned object
#' @param output_filename name of file to be saved
#'
#' @return saves partition file
#' @importFrom dplyr summarise group_by
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' partition_file(hmm_result,"hmm_partitions")
save_partition_file<-function(hmm_result,output_filename){

  # Get the classification from hmm_result object
  classification<-hmm_result$classification

  #RLE is Run length encoding to count each repeated element in the vector
  run_length=rle(hmm_result$classification)

  #Cummulative sum of runnning lengths to get the site number where each class is beginning from
  site=cumsum(run_length$lengths)

  #Append 0 to starting of the list to include beginning of sites
  site=c(0,site)

  #List of `from` sites
  from_site=head(site,-1)+1

  #List of `to` sites
  to_site=tail(site,-1)

  #grab class numbers in order of occurrence
  classes=as.numeric(gsub("C","",run_length$values))

  #create a dataframe similar to transition table
  partition_table=data.frame(from_site,to_site,classes)

  #create a new column combining `from` site anf `to` site
  partition_table$from_to_site=paste(partition_table$from_site,partition_table$to_site,sep='-')

  #creating tibble to group by each class
  partition_tibble<-partition_table %>%
                    group_by(classes) %>%
                    summarise(lst= paste0(from_to_site, collapse = " "))

  #start writing the file
  sink(output_filename)
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

  print(paste("Input files: ",hmm_result$site_input_file,hmm_result$aln_input_file))
  print(paste("Number of sites: ",nrow(hmm_result$data)))
  print(paste("Number of classes: ",ncol(hmm_result$data)))
  print("Sites in each class: ")
  print(table(hmm_result$classification))
  print(paste("Number of transitions: ",length(rle(hmm_result$classification)$value)-1))
  print(paste("Algorithm used to determine the sequence: ",hmm_result$algorithm))

}

