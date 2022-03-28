#' Predict Class function
#'
#' @param sitein relative path of the input .sitelh/.siteprob file
#' @param alninfo relative path of the alignment file
#' @param model Choose model 1,2,3 or 4. Default=4
#' @param iter maximum iterations. Default=10000
#'
#' @return viterbi path, posterior path and convergence truth as a list
#' @importFrom aphid train Viterbi posterior
#' @importFrom testit has_warning
#' @export
#'
#' @examples
#' pred=predict_class("./iqtree/sample.siteprob","./iqtree/sample.alninfo",3)
#' v.pred<- pred[[1]];p.pred<-pred[[2]];conv<-pred[[3]]

predict_class <- function(sitein,alninfo,model,iter){
  if(missing(model)) {
    model=4
  }
  if(missing(iter)) {
    iter=10000
  }

  tab=read.table(alninfo,header=TRUE)
  data=read.table(sitein,header=FALSE,fill=TRUE)

  if(data[1,2]=="LnL"){
    numTrees=(ncol(data)-2)
    colnames(data)<-c("site","LnL",paste("LnLW_",1:numTrees,sep=''))
    data=data[-1,]
    rownames(data)<-NULL
  } else if(data[1,2]=="p1"){
    numTrees=(ncol(data)-1)
    colnames(data) <- data[1,]
    data=data[-1,]
  } else{
    print("Invalid site info file")
  }

  states <- c("Begin",paste("C",1:numTrees,sep=''))

  seq=paste("E",1:numTrees,sep='')[apply(data[,(ncol(data)-numTrees+1):ncol(data)], 1, which.max)]


  # Initial HMM
  ### Define the transition probability matrix
  A <- matrix(c(rep(.01/(numTrees-1),(numTrees+1)**2)),nrow=numTrees+1,ncol=numTrees+1,byrow = TRUE)
  diag(A) <- .99
  A[1,]<-1/numTrees
  A[,1]<-0
  dimnames(A) <- list(from = states, to = states)

  ### Define the emission probability matrix
  if (model==1){
    residues <- paste("E",1:(numTrees),sep='')
    E=matrix(c(rep(.5/(numTrees-1),numTrees)),nrow=numTrees,ncol=numTrees,byrow = TRUE)
  } else if(model==2){
    residues <- paste("E",1:(numTrees+1),sep='')
    seq[tab$Stat=='C']=paste("E",numTrees+1,sep='')
    E=matrix(c(rep(.5/(numTrees),numTrees)),nrow=numTrees,ncol=numTrees+1,byrow = TRUE)
  } else if(model==3){
    residues <- paste("E",1:(numTrees+2),sep='')
    seq[tab$Stat=='U']=paste("E",numTrees+2,sep='')
    seq[tab$Stat=='C']=paste("E",numTrees+1,sep='')
    E=matrix(c(rep(.5/(numTrees+1),numTrees)),nrow=numTrees,ncol=numTrees+2,byrow = TRUE)
  } else if(model==4){
    residues <- paste("E",1:(numTrees+3),sep='')
    seq[tab$Stat=='S']=paste("E",numTrees+3,sep='')
    seq[tab$Stat=='U']=paste("E",numTrees+2,sep='')
    seq[tab$Stat=='C']=paste("E",numTrees+1,sep='')
    E=matrix(c(rep(.5/(numTrees+2),numTrees)),nrow=numTrees,ncol=numTrees+3,byrow = TRUE)
  } else{
    print("Invalid model selected")
  }
  diag(E) <- .5
  dimnames(E) <- list(states = states[-1], residues = residues)
  print(table(seq))

  ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")
  conv<-TRUE

  # Baum-Welch
  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=iter,logspace = FALSE,cores="autodetect",quiet=TRUE))
  if (warn==TRUE){
    print("Not converged")
    conv<-FALSE
    #bw <- train(hmm,seq,method = "Viterbi",logspace = FALSE,quiet=TRUE)
  }

  # get viterbi path
  viterbi = Viterbi(bw,seq)
  viterbi.path=rownames(bw$E)[viterbi$path + 1]

  #get posterior highest state
  post.prob = posterior(bw,seq)
  post.path=tail(states,numTrees)[apply(post.prob, 2, which.max)]
  return(list(viterbi.path,post.path,conv,data[,(ncol(data)-numTrees+1):ncol(data)],bw))
}

#' Predict Class Mixed function
#'
#' @param sitein relative path of the input .sitelh/.siteprob file
#' @param alninfo relative path of the alignment file
#' @param iter maximum iterations. Default=10000
#'
#' @return viterbi path, posterior path and convergence truth as a list
#' @importFrom aphid train Viterbi posterior
#' @importFrom testit has_warning
#' @export
#'
#' @examples
#' pred=predict_class("data100.2t.10k.fa.siteprob","data100.2t.10k.fa.siteprob",3)
#' v.pred<- pred[[1]];p.pred<-pred[[2]];conv<-pred[[3]]
#'
predict_class_mixed <- function(sitein,alninfo,switch,iter){

  if(missing(switch)) {
    switch=1000
  }
  if(missing(iter)) {
    iter=10000
  }

  tab=read.table(alninfo,header=TRUE)
  data=read.table(sitein,header=FALSE,fill=TRUE)

  if(data[1,2]=="LnL"){
    numTrees=(ncol(data)-2)
    colnames(data)<-c("site","LnL",paste("LnLW_",1:numTrees,sep=''))
    data=data[-1,]
    rownames(data)<-NULL
  } else if(data[1,2]=="p1"){
    numTrees=(ncol(data)-1)
    colnames(data) <- data[1,]
    data=data[-1,]
  } else{
    print("Invalid site info file")
  }

  states <- c("Begin",paste("C",1:numTrees,sep=''))

  seq=paste("E",1:numTrees,sep='')[apply(data[,(ncol(data)-numTrees+1):ncol(data)], 1, which.max)]

  residues <- paste("E",1:(numTrees),sep='')


  # Initial HMM
  ### Define the transition probability matrix
  A <- matrix(c(rep(.01/(numTrees-1),(numTrees+1)**2)),nrow=numTrees+1,ncol=numTrees+1,byrow = TRUE)
  diag(A) <- .99
  A[1,]<-1/numTrees
  A[,1]<-0
  dimnames(A) <- list(from = states, to = states)
  ### Define the emission probability matrix
  E=matrix(c(rep(.5/(numTrees-1),numTrees)),nrow=numTrees,ncol=numTrees,byrow = TRUE)
  diag(E) <- .5
  dimnames(E) <- list(states = states[-1], residues = residues)
  ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")
  conv<-TRUE

  # Baum-Welch
  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=switch,logspace = FALSE,cores="autodetect",quiet=TRUE))
  if (warn==TRUE){
    print("model2")
    residues <- paste("E",1:(numTrees+1),sep='')
    seq[tab$Stat=='C']=paste("E",numTrees+1,sep='')
    E=matrix(c(rep(.5/(numTrees),numTrees)),nrow=numTrees,ncol=numTrees+1,byrow = TRUE)
    diag(E) <- .5
    dimnames(E) <- list(states = states[-1], residues = residues)
    hmm <- structure(list(A = A, E = E), class = "HMM")
    #print(hmm$A)
    #print(hmm$E)
    warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=switch,logspace = FALSE,cores="autodetect",quiet=TRUE))
  }
  if (warn==TRUE){
    print("model3")
    residues <- paste("E",1:(numTrees+2),sep='')
    seq[tab$Stat=='U']=paste("E",numTrees+2,sep='')
    seq[tab$Stat=='C']=paste("E",numTrees+1,sep='')
    E=matrix(c(rep(.5/(numTrees+1),numTrees)),nrow=numTrees,ncol=numTrees+2,byrow = TRUE)
    diag(E) <- .5
    dimnames(E) <- list(states = states[-1], residues = residues)
    hmm <- structure(list(A = A, E = E), class = "HMM")
    #print(hmm$A)
    #print(hmm$E)
    warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=switch,logspace = FALSE,cores="autodetect",quiet=TRUE))
  }
  if (warn==TRUE){
    print("model4")
    residues <- paste("E",1:(numTrees+3),sep='')
    seq[tab$Stat=='S']=paste("E",numTrees+3,sep='')
    seq[tab$Stat=='U']=paste("E",numTrees+2,sep='')
    seq[tab$Stat=='C']=paste("E",numTrees+1,sep='')
    E=matrix(c(rep(.5/(numTrees+2),numTrees)),nrow=numTrees,ncol=numTrees+3,byrow = TRUE)
    diag(E) <- .5
    dimnames(E) <- list(states = states[-1], residues = residues)
    hmm <- structure(list(A = A, E = E), class = "HMM")
    warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=iter,logspace = FALSE,cores="autodetect",quiet=TRUE))
  }
  if (warn==TRUE){
    print("Not converged")
    conv<-FALSE
    #bw <- train(hmm,seq,method = "Viterbi",logspace = FALSE,quiet=TRUE)
  }

  # get viterbi path
  viterbi = Viterbi(bw,seq)
  viterbi.path=rownames(bw$E)[viterbi$path + 1]

  #get posterior highest state
  post.prob = posterior(bw,seq)
  post.path=tail(states,numTrees)[apply(post.prob, 2, which.max)]
  return(list(viterbi.path,post.path,conv,data[,(ncol(data)-numTrees+1):ncol(data)],bw))
}

