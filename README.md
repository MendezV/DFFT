# Density-Funtional Fluctuation Theory MLE  (DFFT_MLE)


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
+ [MLE](#mle-demo)
- [Instructions for Use](#instructions-for-use)
+ [MLE](#mle-use)
- [License](./LICENSE)
- [Issues](https://github.com/neurodata/mgc/issues)
- [Pseudocode](#pseudocode)
- [Matlab code](https://github.com/neurodata/mgc-paper/tree/master/Code/MGC)
- [Citation](#citation)
# Overview

The code in this repository can be used to analyze the collective motion of stationary crowds. In particular, input to the code is a csv file in which rows correspond to time frames and columns correspond to the number of individuals inside each bin (quadrats) o

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.


# System Requirements

## Hardware Requirements

The `MGC` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 25 Mbps.

## Software Requirements

### OS Requirements

This package is supported for *Linux* operating systems. The package has been tested on the following systems:

Linux: Ubuntu 16.04  
Mac OSX:  
Windows:  

Before setting up the `MGC` package, users should have `R` version 3.4.0 or higher, and several packages set up from CRAN.

#### Installing R version 3.4.2 on Ubuntu 16.04

the latest version of R can be installed by adding the latest repository to `apt`:

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

which should install in about 20 seconds.

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

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/neurodata/mgc/issues). 

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

### Help

Instructions for help using the `MGC` statistic can be invoked from the `R` terminal window:

```
help("mgc.sample")
```

which produces the `man` page:

```
MGC Sample

Description

The main function that computes the MGC measure between two datasets: It first computes all local correlations, then use the maximal statistic among all local correlations based on thresholding.

Usage

mgc.sample(A, B, option = "mgc")
Arguments

+ A	is interpreted as:
a [nxn] distance matrix, A is a square matrix with zeros on diagonal
a [nxd] data matrix, Otherwise
+ B	is interpreted as:
a [nxn] distance matrix, B is a square matrix with zeros on diagonal
a [nxd] data matrix, Otherwise
+ option='mgc'	is a string that specifies which global correlation to build up-on.
+ 'mgc' use the MGC global correlation.
+ 'dcor' use the dcor global correlation.
+ 'mantel' use the mantel global correlation.
+ 'rank' use the rank global correlation.

and Returns:
+ statMGC is the sample MGC statistic within [-1,1];

+ localCorr consists of all local correlations by double matrix index;

+ optimalScale the estimated optimal scale in matrix single index.

Author(s)

C. Shen
```

# Pseudocode

Pseudocode for the methods employed in the `mgc` package can be found on the [arXiv - MGC](https://arxiv.org/abs/1609.05148) in Appendix C (starting on page 30).

# Citation

For citing code or the paper, please use the citations found in [citation.bib](./citation.bib).
