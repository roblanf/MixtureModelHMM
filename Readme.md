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

## Google Colab notebook

If you have your files, and just want the answer without using R, you can use the google colab notebook here:



## Worked example

### Input data

Let's start with an alignment of 32 species and 10,000 sites. This is the `example.phy` alignment in this repository, which you can download at [this link](https://raw.githubusercontent.com/roblanf/MixtureModelHMM/master/example.phy)

### Getting the right files from an IQ-TREE mixture model analysis

The MixtureModelHMM package will work with the output of _any_ mixture model from IQ-TREE. For this example, we'll focus on the [GHOST model](https://academic.oup.com/sysbio/article/69/2/249/5541793), which is a model that fits mixtures of branch lengths to phylogenetic datasets. To keep this worked example simple, I simulated the `example.phy` dataset above under a GHOST model with 4 classes, and a GTR substitution model. 

> NB: when you are running mixture models on your own data (be that GHOST models or any other mixture model) you should pay close attention to selecting a good model. I used at `GTR+H4` model here because I _simulated_ the data under that model. Usually you wouldn't know the best model, and would want to spend some time figuring out the best avaialable model from your data. I keep it simple here so that we can focus on running the HMM.

The `MixtureModelHMM` package requires a couple of additional files from IQ-TREE that provide information on the likelihood of each site under each model class (the `.sitelh` file) and on some additional features of each site (the `.alninfo` file). So, to use the `MixtureModelHMM` package you'll need to add a couple of things to your usual IQ-TREE command line (this analysis might take a few minutes):

```
iqtree -s example.phy -m GTR+H4 -wslm -alninfo
```

The options we used above are:

* `-s` points to the input alignment file
* `-m` specifies the model. Here we've chosen a GTR model with a 4-class GHOST model (`+H4`), because I know that's the best model
* `-wslm` writes out the site likelihoods to the `.sitelh` file
* `-alninfo` writes out a additional information about sites to the `.alninfo` file, for example which sites are constant, constant but ambiguous, uninformative, or have equal parsimony scores.

Now we have all the files we need for our HMM analysis:

* `example.sitelh`: contains site likelihoods
* `example.alninfo`: contains all the additional information on alignment sites that the HMM needs

### Plotting raw site likelihoods with `plot_scatter()`

The point of the HMM is to take site likelihoods, and leverage the fact that neighbouring sites in the alignment are likely to belong to the same class.

Before running the HMM, it's a really good idea to just look at the raw site likelihoods from IQ-TREE, to see for yourself if neighbouring sites really do have similar class assignments. For this you can use the `plot_scatter()` function as follows:

```{r}
plot_scatter("example.sitelh")
```

Which will give you the following plot:



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


