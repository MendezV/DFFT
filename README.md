# Density-Funtional Fluctuation Theory MLE  (DFFT_MLE)


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for Use](#instructions-for-use)

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

Linux: Ubuntu 16.04  
Windows:  --

#### Other Requirements

Prior to installation and implementation of the code in this repository, one should have one of the following versions of either MATLAB or GNU Octave installed. 

-MATLAB R2017b (9.3.0.713579), 64-bit together with Java 1.8.0_121-b13 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

or,

-GNU Octave, version 4.2.1.

Our code was developed primarily using the former.  

#### Toolboxes and versions on MATLAB

The complete list of toolboxes and versions that where installed at the time of developing the code is:
```
MATLAB                                                Version 9.3         (R2017b)
Simulink                                              Version 9.0         (R2017b)
Aerospace Blockset                                    Version 3.20        (R2017b)
Aerospace Toolbox                                     Version 2.20        (R2017b)
Antenna Toolbox                                       Version 3.0         (R2017b)
Audio System Toolbox                                  Version 1.3         (R2017b)
Automated Driving System Toolbox                      Version 1.1         (R2017b)
Bioinformatics Toolbox                                Version 4.9         (R2017b)
Communications System Toolbox                         Version 6.5         (R2017b)
Computer Vision System Toolbox                        Version 8.0         (R2017b)
Control System Toolbox                                Version 10.3        (R2017b)
Curve Fitting Toolbox                                 Version 3.5.6       (R2017b)
DO Qualification Kit                                  Version 3.4         (R2017b)
DSP System Toolbox                                    Version 9.5         (R2017b)
Database Toolbox                                      Version 8.0         (R2017b)
Datafeed Toolbox                                      Version 5.6         (R2017b)
Econometrics Toolbox                                  Version 4.1         (R2017b)
Embedded Coder                                        Version 6.13        (R2017b)
Filter Design HDL Coder                               Version 3.1.2       (R2017b)
Financial Instruments Toolbox                         Version 2.6         (R2017b)
Financial Toolbox                                     Version 5.10        (R2017b)
Fixed-Point Designer                                  Version 6.0         (R2017b)
Fuzzy Logic Toolbox                                   Version 2.3         (R2017b)
Global Optimization Toolbox                           Version 3.4.3       (R2017b)
HDL Coder                                             Version 3.11        (R2017b)
IEC Certification Kit                                 Version 3.10        (R2017b)
Image Acquisition Toolbox                             Version 5.3         (R2017b)
Image Processing Toolbox                              Version 10.1        (R2017b)
Instrument Control Toolbox                            Version 3.12        (R2017b)
LTE HDL Toolbox                                       Version 1.0         (R2017b)
LTE System Toolbox                                    Version 2.5         (R2017b)
MATLAB Coder                                          Version 3.4         (R2017b)
MATLAB Compiler                                       Version 6.5         (R2017b)
MATLAB Compiler SDK                                   Version 6.4         (R2017b)
MATLAB Distributed Computing Server                   Version 6.11        (R2017b)
MATLAB Report Generator                               Version 5.3         (R2017b)
Mapping Toolbox                                       Version 4.5.1       (R2017b)
Model Predictive Control Toolbox                      Version 6.0         (R2017b)
Neural Network Toolbox                                Version 11.0        (R2017b)
Optimization Toolbox                                  Version 8.0         (R2017b)
Parallel Computing Toolbox                            Version 6.11        (R2017b)
Partial Differential Equation Toolbox                 Version 2.5         (R2017b)
Phased Array System Toolbox                           Version 3.5         (R2017b)
Polyspace Bug Finder                                  Version 2.4         (R2017b)
Polyspace Code Prover                                 Version 9.8         (R2017b)
Powertrain Blockset                                   Version 1.2         (R2017b)
RF Blockset                                           Version 6.1         (R2017b)
RF Toolbox                                            Version 3.3         (R2017b)
Risk Management Toolbox                               Version 1.2         (R2017b)
Robotics System Toolbox                               Version 1.5         (R2017b)
Robust Control Toolbox                                Version 6.4         (R2017b)
Signal Processing Toolbox                             Version 7.5         (R2017b)
SimBiology                                            Version 5.7         (R2017b)
SimEvents                                             Version 5.3         (R2017b)
Simscape                                              Version 4.3         (R2017b)
Simscape Driveline                                    Version 2.13        (R2017b)
Simscape Electronics                                  Version 2.12        (R2017b)
Simscape Fluids                                       Version 2.3         (R2017b)
Simscape Multibody                                    Version 5.1         (R2017b)
Simscape Power Systems                                Version 6.8         (R2017b)
Simulink 3D Animation                                 Version 7.8         (R2017b)
Simulink Check                                        Version 4.0         (R2017b)
Simulink Coder                                        Version 8.13        (R2017b)
Simulink Control Design                               Version 5.0         (R2017b)
Simulink Coverage                                     Version 4.0         (R2017b)
Simulink Design Optimization                          Version 3.3         (R2017b)
Simulink Design Verifier                              Version 3.4         (R2017b)
Simulink Desktop Real-Time                            Version 5.5         (R2017b)
Simulink Report Generator                             Version 5.3         (R2017b)
Simulink Requirements                                 Version 1.0         (R2017b)
Simulink Test                                         Version 2.3         (R2017b)
Stateflow                                             Version 9.0         (R2017b)
Statistics and Machine Learning Toolbox               Version 11.2        (R2017b)
Symbolic Math Toolbox                                 Version 8.0         (R2017b)
System Identification Toolbox                         Version 9.7         (R2017b)
Text Analytics Toolbox                                Version 1.0         (R2017b)
Tracking and Sensor Fusion Toolbox                    Version 1.0         (R2017b)
Trading Toolbox                                       Version 3.3         (R2017b)
WLAN System Toolbox                                   Version 1.4         (R2017b)
Wavelet Toolbox                                       Version 4.19        (R2017b)
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

In any case, the only requirement to use the code is adding the directories in this repository to your local path whether you are working on MATLAB or Octave. The time it takes to download the required files is tipically less than one second. 


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

-```CovMatPoiss```: a (Nbins) square, positive, symmetric, invertible matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators for the naive model. 

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

-```Vmle```: a  vector with 50 entries that corresponds to the vexation that came out of the MLE (if gauge!=0 the average potential is set to zero)

-```fmleError```: a vector with 8 entries that corresponds to the diagonal of the covariance matrix, that ammounts to the variances for each of the parameters in the f-f sector of the covariance matrix, with zeros for the gauge fixed values

-```VmleError```: a  vector with 50 entries that corresponds to the diagonal of the covariance matrix, that ammounts to the variances for each of the parameters in the V-V sector of the covariance matrix

-```CovMatmle ```(if gauge=0): a 56x56 square, positive, symmetric, invertible Matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators. (if there was no gauge fix)

-```CovMatmle``` (if gauge=1): a 57x57 square, positive, symmetric, invertible Matrix that corresponds to the covariance matrix of the asymptotic gaussian distribution for the ML estimators after performing the gauge transformation which ammounts to performing a similarity transformation to the covriance matrix also we append to the asymptotic covariance the error and covariances of the parameter f(1) that was previously fixed but now has error. (if there was a gauge fix)

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
[predictedAv, errorbars, muDFT]=predav(fmle,VPoiss,CovMatmle,CovMatmle,150,'True',0);
```
The output to MATLAB workspace should be :

-```predictedAv``` a vector of size Nbins (the number of bins in the system) in which each entry is the average number of flies as predicted from the DFFT model implemented with a binary search that determines the chemical potential mu that fixes the total number of flies (vector is sorted in the same order as the occupations matrix first dimension)

-```errorbars```: a vector of size Nbins (the number of bins in the system) in which each entry is the one sigma error on the average number of flies as predicted from the DFFT model for each bin(vector is sorted in the same order as the occupations matrix first dimension)

-```muDFT```: Chemical potential that fixes the average total number of flies in the system, for the DFFT model. 



