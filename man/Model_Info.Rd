\name{Model_Info}
\alias{Model_Info}

\title{Details about the Models}

\description{
This project implements Baum-Welch algorithm(B-W) to predict trees/class from output file of Mixtures models such as MAST,GHOST or Q-matrix model.\cr
The input file for B-W algorithm is .sitelh/.siteprob along with .alninfo output files from iqtree.\cr
Project is created using R, depends on R packages "aphid" and "testit".\cr
For plots R packages such as "tidyverse" and "reshape2" are required.\cr
The output is vector indicating the tree/class for a given site.


Model 1 - The tree with highest probability at a given site is considered, and converted into a sequence for training of B-W algorithm.\cr
Model 2 - The tree with highest probability at a given site and constant sites are considered, and converted into a sequence for training of B-W algorithm.\cr
Model 3 - The tree with highest probability at a given site constant sites and non-informative sites are considered, and converted into a sequence for training of B-W algorithm.\cr
Model 4 - The tree with highest probability at a given site constant sites, non-informative sites and same parsimony sites are considered, and converted into a sequence for training of B-W algorithm.\cr
}

