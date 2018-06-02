# Density-Funtional Fluctuation Theory MLE  (DFFT_MLE)


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for Use](#instructions-for-use)

# Overview

The code in this repository can be used to analyze the collective motion of stationary crowds. The study of these crowds is done by binning the system and counting the number of individuals inside each bin at a given instant in time. In particular, the input to the code in this repository is a MAT-file containing a structure array with a single field. This field is a matrix in which rows correspond to time frames and columns correspond to each bin (labeled by an integer) . As such, the (i,j) entry of this matrix corresponds to the number of individuals inside bin j on frame i. With this matrix, the code returns two arrays that caracterize the density-distributions of the crowd among the bins. First, the code returns f(N) (frustration) a vector that characterizes inter-agent interactions within a crowd as a function of the packing within a bin. And second, the code returns $V_B$ (Vexation ), a vector that characterizes agent interactions with their environment (takes lower values near more preferable locations). 

Also, included in this repository is the code to make predictions of crowd density-distributions under new circumstances. By mixing and matching V's and f's from different experiments, one is able to reconstruct the density distribution of a crowd in an environment where V is known for a couple of individuals, but now, the number of individuals increases significantly to a known value.

# Repo Contents

- [Matlab_analysis](./Matlab_analysis): `Matlab` code.

	-[Predictions](./Matlab_analysis/Predictions): Prediction of density-distribution by mixing and matching Vexations and frustrations. Uses binary search to determine the chemical potential of the system to fix the number of agents.
	
	-[MLE_DFFT](./Matlab_analysis/MLE_DFFT): Maximum Likelihood Estimation of both location-preference and inter-agent interaction parameters using the DFFT model
	
	-[MLE_Poiss](./Matlab_analysis/MLE_Poiss): Maximum Likelihood Estimation of location-preference parameters using a Poisson model
	
- [Trial_data](./Trial_data): Test data


# System Requirements

## Hardware Requirements
The code was developed and tested using a computer with the following specifications:

Processor: 2,5 GHz Intel Core i5
RAM: 8 GB 1867 MHz DDR3

The times reported below are calculated with these specs.

## Software Requirements

### OS Requirements

The package has been tested on the following systems:

Mac OSX:  High Sierra Version 10.13.4 (17E202)
Linux: Ubuntu 16.04  
Windows:  --

#### Other Requirements

Before prior to installation and implementation of the code in this repository, one should have one of the following versions of either MATLAB or GNU Octave installed. 

-MATLAB R2017b (9.3.0.713579), 64-bit (maci64) 

or,

-GNU Octave, version 4.2.1.

Our code was developed primarily using the former.

#### Package dependencies

Users should install the following packages prior to installing `mgc`, from an `R` terminal:

```
install.packages(c('ggplot2', 'reshape2', 'Rmisc', 'devtools', 'testthat', 'knitr', 'rmarkdown', 'latex2exp', 'MASS'))
```

which will install in about 80 seconds on a recommended machine.

#### Package Versions

The `mgc` package functions with all packages in their latest versions as they appear on `CRAN` on October 15, 2017. Users can check [CRAN snapshot](https://mran.microsoft.com/timemachine/) for details. The versions of software are, specifically:
```
ggplot2: 2.2.1
reshape2: 1.4.2
Rmisc: 1.5
devtools: 1.13.3
testthat: 0.2.0
knitr: 1.17
rmarkdown: 1.6
latex2exp: 0.4.0
MASS: 7.3
```



# Installation Guide

From an `R` session, type:

```
require(devtools)
install_github('neurodata/mgc', build_vignettes=TRUE)  # install mgc with the vignettes
require(mgc)  # source the package now that it is set up
vignette("MGC", package="mgc")  # view one of the basic vignettes
```

The package should take approximately 20 seconds to install with vignettes on a recommended computer. 

# Demo

## MGC Demo

The `mgc` demo can be run as follows:

```
require(mgc)
set.seed(12345)
mgc.sample(mgc::test_xunif, mgc::test_ylin)$statMGC  # test with a linear relationship between x and y
```

the x data provided is by sampling 100 times from a uniform distribution on the interval [-1, 1], and the y data is formed by adding normally distributed error with variance 0.2 (indicating a linear relationship).

and is expected  to produce the following result exactly approximately *instantaneously*:

```
0.891153
```

a more interactive demo is provided in the package vignette:

```
library(mgc)
vignette("MGC", package="mgc")
```




# Instructions for Use

## MGC Use

In the below tutorial, we show the result of `MGC` to determine the relationship between the first (sepal length) and third (petal length) dimensions of the `iris` dataset, which should run in about 2 seconds:

```
library(mgc)
set.seed(12345)
res <- mgc.sample(iris[,1], iris[,3])
mgc.plot.plot_matrix(res$localCorr, title="MGC Corr Map, Sepal Length and Petal Length",
xlab="Sepal Length Neighbors", ylab="Petal Length Neighbors", legend.name = "statMGC")
print(res$statMGC)
```

![image](https://user-images.githubusercontent.com/8883547/32355967-7de64590-c008-11e7-9c3b-e24470fdbdaa.png)

with the following statistic:

```
0.7337225
```

viewing the corr map above we see that the relationship betweel Sepal and Petal Length is somewhat linear.



