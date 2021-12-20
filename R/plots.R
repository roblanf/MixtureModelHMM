#' Pediction Plots
#'
#' @param pred predict_tree/predict_tree_mixed returned object
#'
#' @return input, viterbi and posterior decoding predicted path plot
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes  coord_flip geom_col scale_fill_continuous geom_point geom_smooth
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' plot_mixtures(pred)
plot_mixtrees<-function(pred){

  v.pred<- as.numeric(gsub("T","",pred[[1]]));p.pred<-as.numeric(gsub("T","",pred[[2]]));seq<-apply(pred[[4]], 1, which.max)
  d<-t(data.frame(v.pred,p.pred,seq))
  row.names(d) <- c("v.pred","p.pred","seq")
  d <- data.frame(names = row.names(d), d)
  df2<-d %>%
    melt(id.vars = "names") %>%
    mutate(variable = 1)

  df2 %>%
    ggplot(aes(x = names, y = variable, group = names, fill = value)) +
    scale_fill_continuous(breaks = unique(seq))+
    geom_col() + coord_flip()
}

#' Inital Plots
#'
#' @param pred predict_tree/predict_tree_mixed returned object
#' @param variable post.prob.tree. or log.like.tree. Default= post.prob.tree.
#'
#' @return scatter plot for initial post.prob.tree./log.like.tree.
#' @export
#'
#' @examples
#' plot_scatter(pred,"post.prob.tree.")
plot_scatter<-function(pred,variable){
  if(missing(variable)) {
    variable="post.prob.tree."
  }
  df=pred[[5]]
  numTrees=(ncol(df)-5)/2
  variables=paste(variable,as.character(1:numTrees),sep='')
  dl = df %>% pivot_longer(cols=variables,names_to = "type", values_to = "measure")
  ggplot(dl, aes(x=site, y = measure, colour=type)) + geom_point(size=0.5,alpha=0.1)+geom_smooth(method='loess', span=0.03)
}
