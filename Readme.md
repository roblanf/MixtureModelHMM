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

Let's start with an alignment of 32 species and 10,000 sites. This is the `example.phy` alignment in this repository, which you can download at [this link](https://raw.githubusercontent.com/roblanf/MixtureModelHMM/master/worked_example/example.phy)

### Getting the right files from an IQ-TREE mixture model analysis

The MixtureModelHMM package will work with the output of _any_ mixture model from IQ-TREE. For this example, we'll focus on the [GHOST model](https://academic.oup.com/sysbio/article/69/2/249/5541793), which is a model that fits mixtures of branch lengths to phylogenetic datasets. To keep this worked example simple, I simulated the `example.phy` dataset above under a GHOST model with 4 classes, and a GTR substitution model. 

> NB: when you are running mixture models on your own data (be that GHOST models or any other mixture model) you should pay close attention to selecting a good model. I used at `GTR+H4` model here because I _simulated_ the data under that model. Usually you wouldn't know the best model, and would want to spend some time figuring out the best avaialable model from your data. I keep it simple here so that we can focus on running the HMM.

The `MixtureModelHMM` package requires a couple of additional files from IQ-TREE that provide information on the likelihood of each site under each model class (the `.sitelh` file) and on some additional features of each site (the `.alninfo` file). So, to use the `MixtureModelHMM` package you'll need to add a couple of things to your usual IQ-TREE command line:

```
iqtree -s example.phy -m GTR+H4 -wslm -wspm -alninfo
```

This analysis will take around 10 minutes (depending on your computer). If you'd like to skip ahead, you can download the output files from this analysis from the [worked_example folder](https://github.com/roblanf/MixtureModelHMM/master/worked_example).

The options we used above are:

* `-s` points to the input alignment file
* `-m` specifies the model. Here we've chosen a GTR model with a 4-class GHOST model (`+H4`), because I know that's the best model
* `-wslm` writes out the site likelihoods under each class to the `.sitelh` file
* `-wspm` writes out the site probabilities under each class to the `.siteprob` file
* `-alninfo` writes out a additional information about sites to the `.alninfo` file, for example which sites are constant, constant but ambiguous, uninformative, or have equal parsimony scores.

Now we have all the files we need for our HMM analysis:

* `example.sitelh`: contains site likelihoods
* `example.siteprob`: contains site probabilities
* `example.alninfo`: contains all the additional information on alignment sites that the HMM needs

### Plotting the raw data with `plot_scatter()`

The point of the HMM is to take site likelihoods or probabilities, and leverage the fact that neighbouring sites in the alignment are likely to belong to the same class to determine tracts of an alignment that belong to the same class.

Before running the HMM, it's a really good idea to just look at the raw data from IQ-TREE, to see for yourself if neighbouring sites really do have similar class assignments. For this you can use the `plot_scatter()` function as follows:

```{r}
library("MixtureModelHMM")
plot_scatter("example.phy.siteprob")
```

> NB: you can also plot the site likelihoods in the same way, but you'll see if you try it that these are far less informative!

The command above will give you the following plot:

![a scatter plot showing a clear pattern](https://github.com/roblanf/MixtureModelHMM/blob/master/img/scatter_plot.png)

Each dot in this plot is a single site in the alignment. The y-axis shows the posterior probability that a site belongs to a given class. This plot shows a few things very clearly. First, neighbouring sites really do seem to share class assignments. For example, sites 1 to about 1250 seem to belong to class `p3` (the blue line). Although the data are very noisy. Second, there are a small number of clear transitions - if you look at which class has the highest probability in the smoothed lines, you'll see that it transitions from p3 to p2, then back to p3, then to p1, and so on. 

This is sufficient to show that it's sensible to run an HMM on these data, so let's go ahead and run it.

### Running the HMM with `run_HMM`

Now we know it's sensible to run an HMM, we can run it in R like this:

```{r}
hmm_result <- run_HMM(site_info = "example.phy.sitelh", aln_info = "example.phy.alninfo")
```

the `run_HMM` function takes two files as input:

* `site_info`: here you can pass either the `.sitelh` or the `.siteprob` file from your IQ-TREE analysis. In amost all of our tests and simulations, we found that the HMM performs better using site likelihoods (`.sitelh` file) rather than site posterior probabilities (the `.siteprob` file), so we recommend using it here.
* `aln_info`: here you have to point to a `.alninfo` file from IQ-TREE.

The HMM trys to figure out how to assign every site in the alignment to a class, using the fact that neighbouring tend to be in the same class to help it. One way to think about this is that the HMM is a way of cleaning up the noisy signal in the plot above. Along the way the HMM accounts for issues associated with phylogenetic data, such as the fact that constant sites don't contain much useful information.

Once the HMM has finished, you can see the results like this:

```
plot_predictions(hmm_result)
```

and write a report on the HMM with a lot more information like this:

```
save_report(hmm_result = hmm_result, output_filename = "hmm_report.txt")
```

and we can save the site classifications of the HMM like this:

```
save_file(hmm_result = hmm_result, output_filename = "hmm_output.gzip")
```

## Description of the approach


