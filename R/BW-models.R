#' Predict Class function
#'
#' @param site_info relative path of the input .sitelh/.siteprob file
#' @param aln_info relative path of the alignment file
#' @param model Choose model 1,2,3 or 4. Default=4
#' @param iter maximum iterations. Default=10000
#' @param algorithm to be used to estimate final sequence of classes
#'
#' @return viterbi path, posterior path and convergence truth as a list
#' @importFrom aphid train Viterbi posterior
#' @importFrom testit has_warning
#' @export
#'
#' @examples
#' hmm_result=run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo",model = 3)
#' viterbi_path<- pred[[1]];p.pred<-pred[[2]];conv<-pred[[3]]

<<<<<<< HEAD
run_HMM <- function(site_info,aln_info,model,iter,algorithm){
  if(missing(model)) {
    model=4
  }
  if(missing(iter)) {
    iter=10000
  }
=======
run_HMM <- function(site_info,aln_info,model=4,iter=10000){
>>>>>>> 41b097b5d5bf10845689e4f28b555d90044dba3b

  tab=read.table(aln_info,header=TRUE)
  data=read.table(site_info,header=FALSE,fill=TRUE)

  # add missing names to IQ-TREE output file
  if(data[1,2]=="LnL"){
    numClasses=(ncol(data)-2)
    colnames(data)<-c("site","LnL",paste("LnLW_",1:numClasses,sep=''))
    data=data[-1,]
    rownames(data)<-NULL
  } else if(data[1,2]=="p1"){
    numClasses=(ncol(data)-1)
    colnames(data) <- data[1,]
    data=data[-1,]
  } else{
    print("Invalid site info file")
  }

  states <- c("Begin",paste("C",1:numClasses,sep=''))

  seq=paste("E",1:numClasses,sep='')[apply(data[,(ncol(data)-numClasses+1):ncol(data)], 1, which.max)]


  # Initial HMM
  ### Define the transition probability matrix
  A <- matrix(c(rep(.01/(numClasses-1),(numClasses+1)**2)),nrow=numClasses+1,ncol=numClasses+1,byrow = TRUE)
  diag(A) <- .99
  A[1,]<-1/numClasses
  A[,1]<-0
  dimnames(A) <- list(from = states, to = states)

  ### Define the emission probability matrix
  if (model==1){
    print("Model 1 - same number of class and emissions")
    residues <- paste("E",1:(numClasses),sep='')
    E=matrix(c(rep(.5/(numClasses-1),numClasses)),nrow=numClasses,ncol=numClasses,byrow = TRUE)
  } else if(model==2){
    print("Model 2 - 1 additional emission for constant site")
    print("The probability of this addtional emission is equal for all classes")
    residues <- paste("E",1:(numClasses+1),sep='')
    seq[tab$Stat=='C']=paste("E",numClasses+1,sep='')
    E=matrix(c(rep(.5/(numClasses),numClasses)),nrow=numClasses,ncol=numClasses+1,byrow = TRUE)
  } else if(model==3){
    print("Model 3 - 2 additional emission for constant site and non-informative sites")
    print("The probability of these addtional emission is equal for all classes")
    residues <- paste("E",1:(numClasses+2),sep='')
    seq[tab$Stat=='U']=paste("E",numClasses+2,sep='')
    seq[tab$Stat=='C']=paste("E",numClasses+1,sep='')
    E=matrix(c(rep(.5/(numClasses+1),numClasses)),nrow=numClasses,ncol=numClasses+2,byrow = TRUE)
  } else if(model==4){
    print("Model 4 - 3 additional emission for constant site, non-informative sites and same parsimony")
    print("The probability of these addtional emission is equal for all classes")
    residues <- paste("E",1:(numClasses+3),sep='')
    seq[tab$Stat=='S']=paste("E",numClasses+3,sep='')
    seq[tab$Stat=='U']=paste("E",numClasses+2,sep='')
    seq[tab$Stat=='C']=paste("E",numClasses+1,sep='')
    E=matrix(c(rep(.5/(numClasses+2),numClasses)),nrow=numClasses,ncol=numClasses+3,byrow = TRUE)
  } else{
    print("Invalid model selected")
  }
  diag(E) <- .5
  dimnames(E) <- list(states = states[-1], residues = residues)
  print("Sites emitted from each class")
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
  post.path=tail(states,numClasses)[apply(post.prob, 2, which.max)]
  return(list(viterbi.path,post.path,conv,data[,(ncol(data)-numClasses+1):ncol(data)],bw))
}
