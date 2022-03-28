## Table of contents
* [General info](#general-info)
* [Models](#models)
* [Installation and Usage](#Installation)
* [Tutorial](#Tutorial)

## General info
This project implements Baum-Welch algorithm to predict trees/class from output file of Mixtures models such as MAST, GHOST or Q-matrix model.\
The input file for B-W algorithm is .sitelh/.siteprob along with .alninfo output files from iqtree.\
Project is created using R, depends on R packages "aphid", "tidyverse", "reshape2" and "testit".\
The output is vector indicating the tree/class for a given site.\
Supports initial scatter plots and prediction plots.

## Models
Model 1 - The tree with highest probability at a given site is considered, and converted into a sequence for training of B-W algorithm.\
Model 2 - The tree with highest probability at a given site and constant sites are considered, and converted into a sequence for training of B-W algorithm.\
Model 3 - The tree with highest probability at a given site, constant sites and non-informative sites are considered, and converted into a sequence for training of B-W algorithm.\
Model 4 - The tree with highest probability at a given site, constant sites, non-informative sites and same parsimony sites are considered, and converted into a sequence for training of B-W algorithm.\
Mix Model - Starts from model 1 and moves to next model if the current B-W model doesn't converge in 100 iterations.

## Installation and Usage
```
devtools::install_github("rahilvora9/PredMixtrees")
library("PredMixtrees")
help(package="PredMixtrees")
```

## Tutorial

Using (IQ-TREE)[https://github.com/iqtree/iqtree2] to generate output files from NEXUS partition input file. The output files are stored in IQTREE_Outfiles folder under the name 'sample_file'. According to the type of file 'sample_files' is given appropriate extension.
```iqtree -s data/sample.nex -m GTR+F+H4 -pre ./IQTREE_Outfiles/sample_file -wspm -wslm -alninfo -redo -nt AUTO```

Creating scatter plots to generate initial log-likelihood/posterior probabilities for each site.
Input file could be either a .sitelh or .siteprob depending which parameter you would like to use for boundary prediction of classes.
```plot_scatter("./IQTREE_Outfiles/sample_file.sitelh")```

To predict class boundaries:
Input file is alignment information(.alninfo) and either a .sitelh or .siteprob depending which parameter you would like to use for boundary prediction of classes. The are other optional arguments which can be specified such as [Models](#models) selection and maximum iterations.

```pred=predict_class("./IQTREE_Outfiles/sample_file.sitelh","./IQTREE_Outfiles/sample_file.alninfo",model=3)
v.pred<- pred[[1]];p.pred<-pred[[2]];conv<-pred[[3]]```

There is other function called predict_class_mixed which starts from model 1 and moves to next model if the current B-W model doesn't converge in user specified iterations.

```pred=predict_class_mixed("./IQTREE_Outfiles/sample_file.sitelh","./IQTREE_Outfiles/sample_file.alninfo", switch=100)
v.pred<- pred[[1]];p.pred<-pred[[2]];conv<-pred[[3]]```

The default value for maximum iterations for both above functions is 10000, while default iterations the `predict_class_mixed` function would move to next model is 1000 iterations. Default model that `predict_class` would use is model 4.

The output is a list where the first element is class prediction of sites done by viterbi algorithm, the second element is class prediction of sites is done using posterior decoding algorithm. The v.pred and p.pred variables are vector associating site with a class number(represented by C, eg C1) chronologically.
The last element is boolean variable to determine whether the algorithm converged or not.

Plotting predictions:
The prediction plots provides sitewise comparison of input and output class association.
```plot_predictions(pred)```

Saving output and report:

Save output vector of class predictions in gz file.
First argument is returned list from prediction function, second is to specify the algorithm such as 'viterbi' or 'posterior' and the final argument is the filename.

```save_file(pred,"viterbi","output")
save_file(pred,"posterior","output")```

Load it back to R:
```v.pred<-system("gzcat viterbi.gz",intern=TRUE)```

Report of final transition and emissions probabilities including proportion of predicted sites for each class can be generated using:

```save_report(pred,"report")```
First argument is returned list from prediction function and second argument is the filename. The file is saves as .txt file
