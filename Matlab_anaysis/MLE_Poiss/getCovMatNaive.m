function [stderrorsNaive,CovMatNaive]=getCovMatNaive(VNaive, MaxPop, Nbins, Tframes)

%%IN
%%-VNaive: a Nbins(number of bins in the system) sized vector that
%%corresponds to the vexation that came out of the MLE for the naive model
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames


%%Calculates the asymptotic Covariance Matrix as the inverse of the fisher
%%information matrix for the naive log-likelihood function

%%OUT
%%-stderrorsNaive: a (Nbins)x1 vector that corresponds to the diagonal of the covariance matrix, corresponding to the
%%variances for each of the parameters that we are estimating for the niave
%%model
%%-CovMatNaive: a (Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators for the naive model. 


%parameters that will be used in the calculation
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system

%%vexation sector of the covariance matrix
z=sum(exp(-VNaive*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%vector that contains the partition function of each bin, size Nbinsx1
NensAv=sum((ones(Nbins,1)*N').*exp(-VNaive*N')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %% the ensemble average according to our model
NsqensAv=sum((ones(Nbins,1)*((N.^2)')).*exp(-VNaive*N')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %% the ensemble average of N^2 according to our model
HessVV=Tframes*diag(NsqensAv-NensAv.*NensAv);

%%final result
Hess=HessVV; %%hessian, equal to the fisher information matrix
CovMatNaive=inv(Hess); %% asymptotic covariance is the inverse of the fisher information
stderrorsNaive=sqrt(diag(CovMatNaive)); %%the factor of two comes in the taylor expansion for the gaussian approximation


end