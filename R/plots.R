#' Prediction Plots
#'
#' @param pred predict_class/predict_class_mixed returned list
#'
#' @return input, viterbi and posterior decoding predicted path plot
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes  coord_flip geom_col geom_point geom_smooth guides guide_legend
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' plot_predictions(hmm_result)
plot_predictions<-function(pred){

  v.pred<- as.numeric(gsub("C","",pred[[1]]));p.pred<-as.numeric(gsub("C","",pred[[2]]));seq<-apply(pred[[4]], 1, which.max)
  d<-t(data.frame(v.pred,p.pred,seq))
  row.names(d) <- c("v.pred","p.pred","input.seq")
  d <- data.frame(names = row.names(d), d)
  df2<-d %>%
    melt(id.vars = "names") %>%
    mutate(Site = 1)

  df2 %>%
    ggplot(aes(x = names, y = Site, group = names, fill = as.factor(value))) +
    geom_col() + coord_flip() + guides(fill=guide_legend(title='Class'))

}

#' Inital Plots
#'
#' @param file relative path of the input sitelh/siteprob file
#'
#' @return scatter plot for initial log-likelihood/posterior probabilities
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth
#' @export
#'
#' @examples
#' plot_scatter("sample_data.sitelh")
#' plot_scatter("sample_data.siteprob")
plot_scatter<-function(site_info){

  data=read.table(site_info,header=FALSE,fill=TRUE)

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
  dl = data %>% pivot_longer(cols=all_of(vars),names_to = "type", values_to = "measure")
  ggplot(dl, aes(x=Site, y = measure, colour=type)) + geom_point(size=0.5,alpha=0.1)+geom_smooth(method='loess', span=0.03)
}

#' Transition Plots
#'
#' @param hmm_result predict_class/predict_class_mixed returned list
#'
#' @return Transition diagram from final probabilities
#' @export
#'
#' @examples
#' plot_hmm_transitions(hmm_result)
#'
plot_hmm_transitions<-function(hmm_result){
  plot(hmm_result[[4]])
}
