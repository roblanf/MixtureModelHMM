\name{Model_Info}
\alias{Model_Info}

\title{Details about the Models}

\description{
This project implements Baum-Welch algorithm to predict trees from output file of Mixtures of trees(MAST) model.\cr
The input file for B-W algorithm is .sitelh output file from iqtree2.\cr
Project is created using R, depends on R packages "aphid" and "testit".\cr
The output is vector indicating the tree for a given site.


Model 1 - The tree with highest probability at a given site is considered, and converted into a sequence for training of B-W algorithm.\cr
Model 2 - The tree with highest probability at a given site and constant sites are considered, and converted into a sequence for training of B-W algorithm.\cr
Model 3 - The tree with highest probability at a given site constant sites and non-informative sites are considered, and converted into a sequence for training of B-W algorithm.\cr
Model 4 - The tree with highest probability at a given site constant sites, non-informative sites and same parsimony sites are considered, and converted into a sequence for training of B-W algorithm.\cr
Mix Model - Starts from model 1 and moves to next model if the current B-W model doesn't converge in 100 iterations.
}
