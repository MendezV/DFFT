function [VPoiss,CovMatPoiss,VPoissError]=extract_params_Poiss(str,plots)
%%'DFFT/Trial_Data/occ_V.mat'
%%%%%%%%%%
%%IN
%%-str: path to a .mat file that unpacks a structure array with a single
%%field. This field corresponds to a matrix of dimensions NbinsxTframes matrix with the number 
%of individuals observed in each bin at each timeframe

%%Extracts the frustration and vexation functions for a single dataset 

%%OUT
%%-VPoiss: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE estimation for the vexation-only
%%model. 
%%-VPoissError: a (Nbins)x1 vector that corresponds to the diagonal of the covariance matrix, corresponding to the
%%variances for each of the vexation parameters, in the vexation-only
%%model. 
%%-CovMatPoiss: a (Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators for the naive model. 

counts=cell2mat(struct2cell(load(str)));
tau=Corr(counts);
counts=counts';

%%
%setting up root values for the gradient search
alpharoot =  0.00010000;
%gradient search with random seed
[VPoiss,CovMatPoiss,VPoissError]=MLEPoiss('random',alpharoot,counts,tau);

%%
%plot params
if plots=='True'

    Nbins=size(counts,1); %total number of bins
    B=((1:(Nbins))-1)'; %vector with bin labels (by integers)

    errorbar(B,VPoiss,VPoissError)
    xlabel('Bin # (B)')
    ylabel('Vexation V(B)')
end
