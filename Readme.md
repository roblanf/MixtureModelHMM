## What is MixtureModelHMM



## Installation
First install the `devtools` package if you don't already have it:

```
install.packages("devtools")
```

Then you can install MixtureModelHMM from this repository like this:

```
library(devtools)
install_github("roblanf/MixtureModelHMM")
```

## QuickStart
If you already have an IQ-TREE analysis where you've produced the following files:

* `mydata.siteprob` or `mydata.sitelh`
* `mydata.alninfo`

You can run the HMM as follows:

```
library("MixtureModelHMM")
hmm_result <- run_HMM(site_info = "mydata.sitelh", aln_info = "mydata.alninfo")
```

Then you can view the key plot like this:

```
plot_predictions(hmm_result)
```

and write a report on the HMM with a lot more information like this:

```
save_file((hmm_result = hmm_result, output_filename = "hmm_report.txt"))
```

## Worked example

Let's start with an alignment of x species and 10,000 sites. You can download the alignment at [this link]()

We suspect for whatever reason that a [GHOST model]() would be a good fit to the data here, so we fit a GHOST model with 4 branch length classes. In addition to the usual IQ-TREE commandline, in this case we need to specify a couple of additional options to make sure we write out the appropriate files for our HMM:

```
iqtree -s alignment.nex -m GTR+I+G+H4 -wslm -alninfo
```

Here's a description options in this IQ-TREE commandline:

* `-s` points to the input alignment file
* `-m` specifies the model. You can use any model of course, but here we've chosen a GTR model, with rates described by invariant sites and a gamma distribution (`+I+G`), and a 4-class GHOST model (`+H4`)
* `-wslm` writes out the site likelihoods, which are the main input variable for our HMM
* `-alninfo` writes out a lot of additional information about the alignment, for example which sites are constant, constant but ambiguous, uninformative, or have equal parsimony scores.

Now we have all the files we need for our HMM analysis, we can move to R to run our HMM.

Before running the HMM, you might want to see what the output of your GHOST model really looks like. Remember this model has 4 classes. That means that the `alignment.sitelh` will contain information on the likelihood of each site under each of the four classes.

```
library("MixtureModelHMM")
hmm_result <- run_HMM(site_info = "mydata.sitelh", aln_info = "mydata.alninfo")
```

Then you can view the key plot like this:

```
plot_predictions(hmm_result)
```

and write a report on the HMM with a lot more information like this:

```
save_file((hmm_result = hmm_result, output_filename = "hmm_report.txt"))
```

## Description of the approach


