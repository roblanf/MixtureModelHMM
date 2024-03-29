#' Run_HMM to Post-processes the output of phylogenetic mixture models using HMMs
#'
#' The function implements the Baum-Welch algorithm for training the HMM and
#' uses dynamic programming to to assign every site in the alignment to a class.
#'
#' @param site_info relative path of the input .sitelh or .siteprob file from your IQ-TREE analysis.
#' In almost all of our tests and simulations, we found that the HMM performs better using site likelihoods(.sitelh file) rather than site posterior probabilities (the .siteprob file), so we recommend using it here.
#' @param aln_info relative path of the alignment file from IQ-TREE analysis.
#' @param model Choose model 1,2,3 or 4. Default=4
#' @param iter The maximum number of EM iterations to carry out before the cycling process is terminated and the partially trained model is returned. Defaults to 100.
#' The maximum change in log likelihood between EM iterations before the cycling procedure is terminated is set to 1E-07.
#' @param algorithm Select dynamic algorithm to predict class boundaries from the trained HMM model parameters.
#'
#' @details Model details
#'
#' Model 1 has no additional emissions. Emissions = classes.
#'
#' Model 2 has one additional emissions for constant sites. Emissions = classes + 1.
#'
#' Model 3 has two additional emissions for constant sites and non-informative sites. Emissions = classes + 2.
#'
#' Model 4 has three additional emissions for constant sites, non-informative sites and same parsimony sites. Emissions = classes + 3.
#'
#' Detail information about approach https://github.com/roblanf/MixtureModelHMM
#'
#' @return an object of class MixtureModelHMM
#' @importFrom aphid train Viterbi posterior
#' @importFrom testit has_warning
#' @export
#'
#' @examples
#' hmm_result = run_HMM(site_info = "mydata.sitelh",aln_info = "mydata.alninfo")

run_HMM <- function(site_info,aln_info,model=4,iter=10000,algorithm="viterbi"){

  print("# Loading input files")
  tab=read.table(aln_info,header=TRUE)
  data=read.table(site_info,header=FALSE,fill=TRUE)

  # add missing names to IQ-TREE output file
  if(data[1,2]=="LnL"){
    #calculate number of classes
    numClasses=(ncol(data)-2)
    colnames(data)<-c("site","LnL",paste("LnLW_",1:numClasses,sep=''))
    data=data[-1,]
    rownames(data)<-NULL
  }
  else if(data[1,2]=="p1"){
    #calculate number of classes
    numClasses=(ncol(data)-1)
    colnames(data) <- data[1,]
    data=data[-1,]
  }
  else{
    print("Invalid site info file")
  }

  #create states
  states <- c("Begin",paste("C",1:numClasses,sep=''))

  #create vector for maximum likelihood or posterior probability
  seq=paste("E",1:numClasses,sep='')[apply(data[,(ncol(data)-numClasses+1):ncol(data)], 1, which.max)]


  # Initial HMM
  print("# Initialising HMM")
  ### Define the transition probability matrix
  A <- matrix(c(rep(.01/(numClasses-1),(numClasses+1)**2)),nrow=numClasses+1,ncol=numClasses+1,byrow = TRUE)
  diag(A) <- .99
  A[1,]<-1/numClasses
  A[,1]<-0
  dimnames(A) <- list(from = states, to = states)

  ### Define the emission probability matrix for each model
  if (model==1){
    print("Model 1 - same number of class and emissions")
    residues <- paste("E",1:(numClasses),sep='')
    E=matrix(c(rep(.5/(numClasses-1),numClasses)),nrow=numClasses,ncol=numClasses,byrow = TRUE)
  }
  else if(model==2){
    print("Model 2 - 1 additional emission for constant site")
    print("The probability of this addtional emissions is equal for all classes")
    residues <- paste("E",1:(numClasses+1),sep='')
    seq[tab$Stat=='C']=paste("E",numClasses+1,sep='')
    E=matrix(c(rep(.5/(numClasses),numClasses)),nrow=numClasses,ncol=numClasses+1,byrow = TRUE)
  }
  else if(model==3){
    print("Model 3 - 2 additional emission for constant site and non-informative sites")
    print("The probability of these addtional emissions is equal for all classes")
    residues <- paste("E",1:(numClasses+2),sep='')
    seq[tab$Stat=='U']=paste("E",numClasses+2,sep='')
    seq[tab$Stat=='C']=paste("E",numClasses+1,sep='')
    E=matrix(c(rep(.5/(numClasses+1),numClasses)),nrow=numClasses,ncol=numClasses+2,byrow = TRUE)
  }
  else if(model==4){
    print("Model 4 - 3 additional emission for constant site, non-informative sites and same parsimony")
    print("The probability of these addtional emissions is equal for all classes")
    residues <- paste("E",1:(numClasses+3),sep='')
    seq[tab$Stat=='S']=paste("E",numClasses+3,sep='')
    seq[tab$Stat=='U']=paste("E",numClasses+2,sep='')
    seq[tab$Stat=='C']=paste("E",numClasses+1,sep='')
    E=matrix(c(rep(.5/(numClasses+2),numClasses)),nrow=numClasses,ncol=numClasses+3,byrow = TRUE)
  }
  else{
    print("Invalid model selected")
  }
  diag(E) <- .5
  dimnames(E) <- list(states = states[-1], residues = residues)
  print("Sites emitted from each class")
  print(table(seq))

  ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")

  # Baum-Welch
  print("# Training HMM")

  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=iter,logspace = FALSE,cores="autodetect",quiet=TRUE))

  # Check if training is converged
  if (warn==TRUE){
    warning("Your HMM analysis did not converge. To address this, you should increase the number of iterations. You used N iterations, and a good starting point given that this didn't converge is to double this to Y iterations.
            You can specify the number of iterations in the run_HMM() function using the `iter` argument", call. = FALSE)
  }

  print(paste("# Running",algorithm,"algorithm"))
  if (algorithm == "viterbi"){
    # get viterbi path
    viterbi = Viterbi(bw,seq)
    classification=rownames(bw$E)[viterbi$path + 1]
  }
  else if(algorithm == "posterior"){
    #get posterior highest state
    post.prob = posterior(bw,seq)
    classification=tail(states,numClasses)[apply(post.prob, 2, which.max)]
  }

  #creating output object
  res<-structure(list(classification=classification,data=data[,(ncol(data)-numClasses+1):ncol(data)],trained_hmm=bw,algorithm=algorithm,
                      site_input_file=site_info,aln_input_file=aln_info),class="MixtureModelHMM")

  res$alignment_plot=plot_predictions(res)

  res$hmm_probabilities=list(class_transition_probabilities=bw$A,class_emission_probabilities=bw$E)

  res$hmm_transition_table=transition_table(res)

  res$initial_scatter_plot=plot_scatter(gsub(".sitelh",".siteprob",site_info))
  return(res)
}
