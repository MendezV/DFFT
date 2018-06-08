# Density-Funtional Fluctuation Theory MLE  (DFFT_MLE)


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Instructions for Use and Demo](#instructions-for-use-and-demo)
- [License](#license)


# Overview

The code in this repository can be used to analyze the collective motion of stationary crowds. The study of these crowds is done by binning the system and counting the number of individuals inside each bin at a given instant in time. In particular, the input to the code in this repository is a MAT-file containing a structure array with a single field. This field is a matrix in which rows correspond to time frames and columns correspond to each bin (labeled by an integer) . As such, the (i,j) entry of this matrix corresponds to the number of individuals inside bin j on frame i. With this matrix, the code returns two arrays that caracterize the density-distributions of the crowd among the bins. First, the code returns f(N) (frustration) a vector that characterizes inter-agent interactions within a crowd as a function of the packing within a bin. And second, the code returns V(B) (Vexation ), a vector that characterizes agent interactions with their environment (takes lower values near more preferable locations). 

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

Mac OSX:  High Sierra Version 10.13.4 Build: 17E202 

Linux: Ubuntu 16.04 LTS  
Windows:  --

#### Other Requirements

Prior to installation and implementation of the code in this repository, one should have one of the following versions of either MATLAB or GNU Octave installed. 

-MATLAB Version: 9.0.0.341360 (R2016a), 64-bit together with Java 1.7.0_60-b19 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

or,

-GNU Octave, version 4.2.1.

Our code was developed primarily using the former.  

#### Toolboxes and versions on MATLAB

The complete list of toolboxes and versions that where installed at the time of developing the code is:
```
MATLAB                                                Version 9.0         (R2016a)
Simulink                                              Version 8.7         (R2016a)
Bioinformatics Toolbox                                Version 4.6         (R2016a)
Control System Toolbox                                Version 10.0        (R2016a)
DSP System Toolbox                                    Version 9.2         (R2016a)
Econometrics Toolbox                                  Version 3.4         (R2016a)
Financial Instruments Toolbox                         Version 2.3         (R2016a)
Financial Toolbox                                     Version 5.7         (R2016a)
Fixed-Point Designer                                  Version 5.2         (R2016a)
Global Optimization Toolbox                           Version 3.4         (R2016a)
Image Processing Toolbox                              Version 9.4         (R2016a)
Neural Network Toolbox                                Version 9.0         (R2016a)
Optimization Toolbox                                  Version 7.4         (R2016a)
Partial Differential Equation Toolbox                 Version 2.2         (R2016a)
Signal Processing Toolbox                             Version 7.2         (R2016a)
Simulink Control Design                               Version 4.3         (R2016a)
Statistics and Machine Learning Toolbox               Version 10.2        (R2016a)
Symbolic Math Toolbox                                 Version 7.0         (R2016a)
Wavelet Toolbox                                       Version 4.16        (R2016a)
```
Note that many of these toolboxes were are not needed for the use of our code. Nevertheles, we include this list for completeness. Also, note that there is no need for any non-standard hardware to use the code in this repository.

# Installation Guide

You can fork this repository to your matlab session using 
```
Source Control -> Manage Files... 
```
And typing 
```
https://github.com/MendezV/DFFT
```
in ```repository path```.

In any case, the only requirement to use the code is adding the directories in this repository to your local path whether you are working on MATLAB or Octave. The time it takes to download the required files is typically less than one second. 


# Instructions for Use and Demo

## DFFT_MLE Demo

The `DFFT_MLE` demo can be run as follows:

```
[fmle,Vmle,CovMatmle,fmleError,VmleError]=extract_params_DFFT('DFFT/Trial_Data/occ_f.mat','True');
```
The data in  ```occ_f.mat``` is a MAT-file containing a structure array with a single field. This field is a matrix in which rows correspond to time frames and columns correspond to each bin (labeled by an integer) . As such, the (i,j) entry of this matrix corresponds to the number of individuals inside bin j on frame i. These entries are measured for a crowd of 135 flies inside a circular chamber with 50 bins. The maximum occupation observed in this data is 7 flies in a bin. The total number of frames was 533.  As such, ```occ_f.mat```  contains a 533x50 matrix. 

The output to MATLAB workspace should be :

-```fmle```: a vector with 8 entries that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)

-```Vmle```: a  vector with 50 entries that corresponds to the vexation that came out of the MLE (if gauge!=0 the average potential is set to zero)

-```fmleError```: a vector with 8 entries that corresponds to the diagonal of the covariance matrix, that ammounts to the variances for each of the parameters in the f-f sector of the covariance matrix, with zeros for the gauge fixed values

-```VmleError```: a  vector with 50 entries that corresponds to the diagonal of the covariance matrix, that ammounts to the variances for each of the parameters in the V-V sector of the covariance matrix

-```CovMatmle ```(if gauge=0): a 56x56 square, positive, symmetric, invertible Matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators. (if there was no gauge fix)

-```CovMatmle``` (if gauge=1): a 57x57 square, positive, symmetric, invertible Matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators after performing the gauge transformation which ammounts to performing a similarity transformation to the covriance matrix also we append to the asymptotic covariance the error and covariances of the parameter f(1) that was previously fixed but now has error. (if there was a gauge fix)

the ```gauge``` parameter is set within the ```extract_params_DFFT.mat``` script, and the default setting is ```gauge=0```.  Additionally, the following plot of fmle and Vmle with their respective error bars shoud appear

![image](https://github.com/MendezV/MLE-for-DFT-master/blob/master/Other/Figures/DEMO.png)
 
 This operation typically lasts less than one second. Moreover, the following output should appear in the command line:
 ```
 Correlation time (in frames)...
 
 tau =
 
 6.5316
```
Which tells time needed to decorrelate the system in units of frames. Also, 
```
Elapsed time is 0.452752 seconds.
```
Tells the time in which  the non-linear conjugate gradients algorithm converged. This time changes with each call of the function because random initial conditions for the search where set by default. However, this is time is generally lower than a second for the data in this repository. Note that this time depends monotonically on the ammount of bins in the system. Finally, the line:

```
counter =

1920
```
tells the number of iterations until the non-linear conjugate gradients algorithm converges. For the data in this repository, the number is typically around 2000-3000 depending on the random initial condition. 
 
## Poiss_MLE Demo

The `Poiss_MLE` demo can be run as follows:

```
[VPoiss,CovMatPoiss,VPoissError]=extract_params_Poiss('DFFT/Trial_Data/occ_V.mat','True');
```
The data in  ```occ_V.mat``` is a MAT-file containing a structure array with a single field. This field is a matrix in which rows correspond to time frames and columns correspond to each bin (labeled by an integer) . As such, the (i,j) entry of this matrix corresponds to the number of individuals inside bin j on frame i. These entries are measured for a single flies inside a square chamber with 49 bins. The total number of frames was 6384.  As such, ```occ_V.mat``` contains a 6384x49 matrix. 

The output to MATLAB workspace should be :

-```VPoiss```:  a  vector with 49 entries  that corresponds to the vexation that came out of the MLE estimation for the vexation-only model. 

-```VPoissError```: a  vector with 49 entries that corresponds to the diagonal of the covariance matrix, corresponding to the variances for each of the vexation parameters, in the vexation-only model. 

-```CovMatPoiss```: a (49x49) square, positive, symmetric, invertible matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators for the naive model. 

Additionally, the following plot of VPoiss with its respective error bars shoud appear

![image](https://github.com/MendezV/MLE-for-DFT-master/blob/master/Other/Figures/DEMO2.png)

This operation typically lasts less than one second. Moreover, the following output should appear in the command line:
```
Correlation time (in frames)...

tau =

4.4359
```
Which tells time needed to decorrelate the system in units of frames. Also, 
```
Elapsed time is 0.143906 seconds.
```
Tells the time in which  the non-linear conjugate gradients algorithm converged. This time changes with each call of the function because random initial conditions for the search where set by default. However, this is time is generally lower than a second for the data in this repository. Note that this time depends monotonically on the ammount of bins in the system. Finally, the line:

```
counter =

1081
```
tells the number of iterations until the non-linear conjugate gradients algorithm converges. For the data in this repository, the number is typically around 800-1500 depending on the random initial condition. 

## Predictions Demo

To test the predictions demo, first run the DFFT_MLE demo with

```
[fmle,Vmle,CovMatmle,fmleError,VmleError]=extract_params_DFFT('DFFT/Trial_Data/occ_exp.mat','True');
```
The data in  ```occ_exp.mat``` is a MAT-file containing a structure array with a single field. This field is a matrix in which rows correspond to time frames and columns correspond to each bin (labeled by an integer) . As such, the (i,j) entry of this matrix corresponds to the number of individuals inside bin j on frame i. These entries are measured for a crowd of 135 flies inside a square chamber with 49 bins. The total number of frames was 15960.  As such, ```occ_exp``` contains a 15960x49 matrix. 

The output to MATLAB workspace should be :

-```fmle```: a vector with 8 entries that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)

-```Vmle```: a  vector with 49 entries that corresponds to the vexation that came out of the MLE (if gauge!=0 the average potential is set to zero)

-```fmleError```: a vector with 8 entries that corresponds to the diagonal of the covariance matrix, that ammounts to the variances for each of the parameters in the f-f sector of the covariance matrix, with zeros for the gauge fixed values

-```VmleError```: a  vector with 49 entries that corresponds to the diagonal of the covariance matrix, that ammounts to the variances for each of the parameters in the V-V sector of the covariance matrix

-```CovMatmle ```(if gauge=0): a 55x55 square, positive, symmetric, invertible Matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators. (if there was no gauge fix)

-```CovMatmle``` (if gauge=1): a 56x56 square, positive, symmetric, invertible Matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators after performing the gauge transformation which ammounts to performing a similarity transformation to the covriance matrix also we append to the asymptotic covariance the error and covariances of the parameter f(1) that was previously fixed but now has error. (if there was a gauge fix)

the ```gauge``` parameter is set within the ```extract_params_DFFT.mat``` script, and the default setting is ```gauge=0```.  Additionally, the following plot of fmle and Vmle with their respective error bars shoud appear

![image](https://github.com/MendezV/MLE-for-DFT-master/blob/master/Other/Figures/DEMO3.png)

This operation typically lasts less than one second. Moreover, the following output should appear in the command line:
```
Correlation time (in frames)...

tau =

42.5114

```
Which tells time needed to decorrelate the system in units of frames. Also, 
```
Elapsed time is 0.185710 seconds.
```
Tells the time in which  the non-linear conjugate gradients algorithm converged. This time changes with each call of the function because random initial conditions for the search where set by default. However, this is time is generally lower than a second for the data in this repository. Note that this time depends monotonically on the ammount of bins in the system. Finally, the line:

```
counter =

687
```
tells the number of iterations until the non-linear conjugate gradients algorithm converges. For the data in this repository, the number is typically around 500-1000 depending on the random initial condition. 

Now that the variables are set, a simple test consists of predicting the average density in each bin for the same ```occ_exp.mat``` data.  This can be done by typing:

```
Nflies=135
gauge=0
[predictedAv, errorbars, muDFT]=predav(fmle,Vmle,CovMatmle,CovMatmle,Nflies,'True',gauge);
```
The output to MATLAB workspace should be :

-```predictedAv``` a vector of size Nbins (the number of bins in the system) in which each entry is the average number of flies as predicted from the DFFT model implemented with a binary search that determines the chemical potential mu that fixes the total number of flies (vector is sorted in the same order as the occupations matrix first dimension)

-```errorbars```: a vector of size Nbins (the number of bins in the system) in which each entry is the one sigma error on the average number of flies as predicted from the DFFT model for each bin(vector is sorted in the same order as the occupations matrix first dimension)

-```muDFT```: Chemical potential that fixes the average total number of flies in the system, for the DFFT model. 

To visualize the average density-distribution of flies within the square chamber one can run the following line 
```
heatmap(reshape(predictedAv,[7,7]));
```
which should output the following plot:

![image](https://github.com/MendezV/MLE-for-DFT-master/blob/master/Other/Figures/DEMO4.png)


