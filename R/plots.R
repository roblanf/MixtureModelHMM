#' Prediction Plots
#'
#' @param hmm_result run_hmm() returned object
#'
#' @return hmm_input, viterbi or posterior decoding predicted path plot
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes  coord_flip geom_col geom_point geom_smooth guides guide_legend
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' plot_predictions(hmm_result)
plot_predictions<-function(hmm_result){

  classification<- as.numeric(gsub("C","",hmm_result$classification));seq<-apply(hmm_result$data, 1, which.max)
  d<-t(data.frame(classification,seq))
  row.names(d) <- c(paste("classification",hmm_result$algorithm),"hmm_input_seq")
  d <- data.frame(names = row.names(d), d)

  df2<-d %>%
    melt(id.vars = "names") %>%
    mutate(Site = 1)

  print("# Making alignment_plot")
  df2 %>%
    ggplot(aes(x = names, y = Site, group = names, fill = as.factor(value))) +
    geom_col() + coord_flip() + guides(fill=guide_legend(title='Class'))

}

#' Inital Plots
#'
#' @param input_filepath file path to a .sitelh or a .siteprob file from IQ-TREE
#' @param span the span of the kernel density plot (default 0.03)
#'
#' @return ggplot object with a scatter plot and smoothed kernel density
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth
#' @export
#'
#' @examples
#' plot_scatter(input_filepath = "sample_data.sitelh")
#' plot_scatter(input_filepath = "sample_data.siteprob")

plot_scatter<-function(input_filepath, span=0.03){

  data=read.table(input_filepath,header=FALSE,fill=TRUE)

  if(data[1,2]=="LnL"){
    numClasses=(ncol(data)-2)
    colnames(data)<-c("Site","LnL",paste("LnLW_",1:numClasses,sep=''))
    data=data[-1,]
    rownames(data)<-NULL
  } else if(data[1,2]=="p1"){
    numClasses=(ncol(data)-1)
    colnames(data) <- data[1,]
    data=data[-1,]
  } else{
    print("Invalid site info file")
  }

  data[] <- lapply(data, function(x) as.numeric(as.character(x)))

  vars=colnames(data[,(ncol(data)-numClasses+1):ncol(data)])

  dl = data %>%
         pivot_longer(cols=all_of(vars),names_to = "type", values_to = "measure")

  ggplot(dl, aes(x=Site, y = measure, colour=type)) +
    geom_point(size=0.5,alpha=0.1) +
    geom_smooth(method='loess', span=span)
}

#' Transition Plots
#'
#' @param hmm_result run_hmm() returned object
#'
#' @return Transition diagram from final probabilities
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' plot_hmm_transitions(hmm_result)
#'
plot_hmm_transitions<-function(hmm_result){
  plot(hmm_result$trained_hmm)
}
