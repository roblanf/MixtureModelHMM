## Table of contents
* [Introduction](#what-is-mixturemodelhmm)
* [Installation](#installation)
* [Quick Start](#quickstart)
* [Example](#worked-example)
* [Approach](#description-of-the-approach)

## What is MixtureModelHMM

MixtureModelHMM implements _Hidden Markov Model_ a probabilistic model used to analyze sequencial data. It is used to post-process output of phylogenetic mixture models from [IQ-TREE](http://www.iqtree.org/). The package implements `Baum-Welch` Algorithm to train the HMM model on given input file. The Baum-Welch algorithm is a dynamic programming approach and a special case of the expectation-maximization algorithm (EM algorithm). Its purpose is to tune the parameters of the HMM, namely the state transition matrix, the emission matrix, and the initial state distribution, such that the model is maximally like the observed data. The input files can be either the site likelihood or site probability file and the alignment information file. Once the HMM is trained the model could be used to determine the final class associated with each sites using dynamic programming algorithm. The default algorithm for predicting class boundries is set to `veterbi`. The output of the function returns an object class consisting of vector for each classes assigned to a given class in a [vector form](#classification) along with a [prediction plot](#alignment_plot) and [hmm transition table](#hmm_transition_table).


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
hmm_result$alignment_plot
```

and write a report on the HMM with a lot more information like this:

```
save_report(hmm_result = hmm_result, output_filename = "hmm_report.Rmd")
```

## Google Colab notebook

If you have your files, and just want the answer without using R, you can use the google colab notebook here:
[MixtureModelHMM_notebook](https://colab.research.google.com/drive/1Kn7nigiNPh25tPqn4mC6e-aDYy4njN60?usp=sharing#scrollTo=4CLAR_gBqw8f)

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

```r
library("MixtureModelHMM")
plot_scatter("example.phy.siteprob")
```

> NB: you can also plot the site likelihoods in the same way, but you'll see if you try it that these are far less informative!

The functions implements `geom_smooth()` with default `span=0.03`.\
`span` can be  adjusted as suitable where smaller numbers produce wigglier lines, larger numbers produce smoother lines.

The command above will give you the following plot:

![a scatter plot showing a clear pattern](https://github.com/roblanf/MixtureModelHMM/blob/master/img/scatter_plot.png)

Each dot in this plot is a single site in the alignment. The y-axis shows the posterior probability that a site belongs to a given class. This plot shows a few things very clearly. First, neighbouring sites really do seem to share class assignments. For example, sites 1 to about 1250 seem to belong to class `p3` (the blue line). Although the data are very noisy. Second, there are a small number of clear transitions - if you look at which class has the highest probability in the smoothed lines, you'll see that it transitions from p3 to p2, then back to p3, then to p1, and so on. 

This is sufficient to show that it's sensible to run an HMM on these data, so let's go ahead and run it.

### Running the HMM with `run_HMM`

Now we know it's sensible to run an HMM, we can run it in R like this:

```r
hmm_result <- run_HMM(site_info = "example.phy.sitelh", aln_info = "example.phy.alninfo")
```

the `run_HMM` function takes two files as input:

* `site_info`: here you can pass either the `.sitelh` or the `.siteprob` file from your IQ-TREE analysis. In almost all of our tests and simulations, we found that the HMM performs better using site likelihoods (`.sitelh` file) rather than site posterior probabilities (the `.siteprob` file), so we recommend using it here.
* `aln_info`: here you have to point to a `.alninfo` file from IQ-TREE.

Other optional parameters:

* `model`: Select a model from 1-4. Defaults to model 4.
* `iter` : The maximum number of EM iterations to carry out before the cycling process is terminated and the partially trained model is returned.       Defaults to 100. The maximum change in log likelihood between EM iterations before the cycling procedure is terminated is set to 1E-07.
* `algorithm`: Select dynamic algorithm to predict class boundaries from the trained HMM model parameters.

The HMM trys to figure out how to assign every site in the alignment to a class, using the fact that neighbouring tend to be in the same class to help it. One way to think about this is that the HMM is a way of cleaning up the noisy signal in the plot above. Along the way the HMM accounts for issues associated with phylogenetic data, such as the fact that constant sites don't contain much useful information.

You can get a quick summary of the hmm result like this:


```r
summary(hmm_result)
```

which in this case will show:

```
[1] "Input files:  example.phy.sitelh example.phy.alninfo"
[1] "Number of sites:  10000"
[1] "Number of classes:  4"
[1] "Sites in each class: "

  C1   C2   C3   C4 
2500 2497 2500 2503 
[1] "Number of transitions:  8"
[1] "Algorithm used to determine the sequence:  viterbi"
```

### The key output from the HMM

Most people will be primarily interested in three things from the HMM. Here they are.

#### `alignment_plot`

```r
hmm_result$alignment_plot
```
Will show you this plot:

![a plot of the hmm output along the alignment](https://github.com/roblanf/MixtureModelHMM/blob/master/img/alignment_plot.png)

This plot shows the alignment sites along the x-axis. On the top panel it shows the maximum likelihood (or posterior probability, if that's the file you used as input) class for each site in your alignment. This is the input data for the HMM. On the bottom panel it shows the output of the HMM - i.e. which classes the HMM has assigned each site to. You can compare this plot to the scatter plot above, and see how the HMM cleans up the noisy signal in the input data.

#### `hmm_transition_table`

The plot above shows that the HMM has a few transitions between classes. If we go from right to left, the classes go: 3, 2, 3, 1, 2, 1, 4, 2, 4. The `hmm_transition_table` shows you exactly which sites the transitions occurred at.

```r
hmm_result$hmm_transition_table
``` 

Will show you this:

```
  site class_from class_to
1 1501         C3       C2
2 2301         C2       C3
3 3301         C3       C1
4 4301         C1       C2
5 5101         C2       C1
6 6601         C1       C4
7 7103         C4       C2
8 8000         C2       C4
```

#### `classification`

`classification` is a vector where you can look up the HMM classification of every single site in your alignment.

For example, if we wanted to see the transition between class 3 and class 2, then we might want to look at the sites around site 1501 (the location of the transition above). We could do that as follows:

```r
hmm_result$classification[1495:1505]
```

This will show you the following:

```
 [1] "C3" "C3" "C3" "C3" "C3" "C3" "C2" "C2" "C2" "C2" "C2"
```

And you can see that at site 1500 the output is `C3`, and it switches to `C2` at site 1501.

### Output options

Obviously you can output anything you like about the HMM from the `hmm_results` object (full details below). 

The `MixtureModelHMM` package has a few functions to write out specific files that we hope are helpful.

#### `save_report()`

This will save a Rmd and pdf file that contains a lot of useful information about your results:

```r
save_report(hmm_result = hmm_result, output_filename = "hmm_report")
```

#### `save_partitioning_scheme()`

If you want to output the information from the HMM as a partitioning scheme in NEXUS format, you can do this:

```r
save_partitioning_scheme(hmm_result = hmm_result, output_filename = "hmm_partitions.nex")
```

### Everything that's in the output

The `hmm_result` object is an object of class 'MixtureModelHMM'. This object contains all the information we think you might ever be interested in.

* `$classification`: a vector the same length as your alignment, which contains the final classification for each site
* `$data`: the input data loaded via the `site_info` argument to `run_HMM()`
* `$trained_hmm`: an object of class `HMM`, which describes the final Hidden Markov Model you trained
* `$algorithm`: the name of the algorithm used to get the final path through the HMM
* `$site_input_file`: the name of the site input file loaded via `site_info` to `run_HMM()`
* `$aln_input_file`: the name of the alignment input file loaded via `aln_info` to `run_HMM()`
* `$alignment_plot`: a plot that shows the input data vs. what the HMM inferred
* `$transition_plot`: a plot that shows the structure of the final HMM
* `$hmm_transition_table`: a table that describes each inferred transition from the HMM


## Description of the approach

#### Flowchart
![Flowchart](https://user-images.githubusercontent.com/11074196/173984167-38b6e15a-abe5-4959-b228-9100c04912c0.png)


We can either use posterior probability for each site(.siteprob file) or log-likelihood for each site(.sitelh file) as an input for `run_HMM()`.

* Step 1: We create a sequence of classes with maximum log-likelihood(or posterior proabilility) for each given site.
          
  $Cx_1,Cx_2,Cx_3,...,Cx_m$\
  each $x_i = {\operatorname{argmax}}\set{Var_1,Var_2,Var_3,...,Var_n}$\
  where $Var_i$ = LnLW_i or p_i from `site_info` file, $n$ = number of classes and $m$ = number of sites\
  The above sequence is used to train Baum-Welch algorithm.

* Step 2: Initialize HMM with initial probabilities to start training.

Hidden Markov Model consists of transition probability and emission probability.\
We define states as classes and emissions as classes along with additional states.
* Initialize transition probabilities 
  * Probability of beginning with any class is equally distributed
  * Probability of a class tranisitioning to same class is 0.99 and the 0.01 is distributed equally to n-1 reaminimg classes.
* Initialize Emission probabilities
  * Probability of a class emitted of itself is 0.5 and the other 0.5 is distibuted among other emissions.
  * The probability of this addtional emissions is equal for all classes.
  * Number of emissions are decided according to model selected:
    * Model 1 has no additional emissions. Emissions = classes
    * Model 2 has one additional emissions for constant sites. Emissions = classes + 1
    * Model 3 has two additional emissions for constant sites and non-informative sites. Emissions = classes + 2
    * Model 4 has three additional emissions for constant sites, non-informative sites and same parsimony sites. Emissions = classes + 3
         
* Step 3: We train the parameters and get final transition/emission probs.
![BW Flowchart](https://user-images.githubusercontent.com/11074196/173843794-029ac6ab-f8b6-45e6-8d67-d4c35a7f26a8.png)

* Step 4: Use the final parameters to predict class boundaries using dynamic programming algorithms - Viterbi or posterior decoding(forward-backward).

