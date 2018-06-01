function [fmle,Vmle,CovMatmle,fmleError,VmleError]=MLE(rootparams,alpharoot,counts,gauge,tau)
%%%%%%%%%%
%%IN
%%-rootparams: a vector of size MaxPop+1 + Nbins (MaxPop+1 values for the frustration,the number of flies 
%%can go from zero to the maximum observed packing in all the bins. And Nbins values for the vexation at each bin )
%%it corresponds to the root of the minimization algorithm 
%%-rootalpha: scalar that corresponds to the root that will be used in
%%the linesearch algorithm that is implemented in the minimization 
%%-counts: a NbinsxTframes matrix with the number of flies observed in each
%%bin at each timeframe
%%-gauge numerical value, if 0 no gauge transformation is applied to the
%%parameters and f(0)=f(1)=0. else the gauge is set so that the average
%%potential is zero and f(1) corresponds to the sum of the previous values
%%for the potential
%%-tau:correlation time corrects for the number of independent time samples


%%finds the minimum of g(V,F)=-logP(V,F|data) by the method of preconditioned conjugate
%%gradients the argument of this function at the minimum gives the maximum
%%likelihood estimates for the parameters in our model. The minimization
%%constraints the values of f(0) and f(1) to be 0 so that the covariance
%%matrix is non-singular and there is a well defined gaussian asymptotic
%%limit for the distribution of estimators

%%OUT
%%-f: a MaxPop+1(MaxPop is the maximum number of flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values (if gauge!=0 the average potential is set to zero and the gauge is adjusted for f(1) with appropiate error propagation)
%%-V: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE (if gauge!=0 the average potential is set to zero)
%%-stderrors: a (MaxPop-1+Nbins)x1 vector that corresponds to the diagonal of the covariance matrix, corresponding to the
%%variances for each of the parameters that we are estimating without the
%%gauge fixed values, which have no uncertainty
%%-CovMat: a (MaxPop-1+Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators. (if there was no gauge fix)
%%-CovMat:a (MaxPop+Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators after performing the gauge transformation
%%which ammounts to performing a similarity transformation to the covriance matrix
%%also we append to the asymptotic covariance the error and covariances of the parameter f(1) that was previously fixed but now has error. (if there was a gauge fix)

more off;

%%%%%%%%
MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
hist=delta';
					  
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system
NexpAv=hist*N; %experimental average number of flies at each bin
Nfac=gamma(ones(Nbins,1)*(N+1)');
					  

%%first iteration is steepest descent with a variable step size 
presgrad=logligrad(rootparams ,  MaxPop , Nbins, Tframes, hist, N, Nfac, NexpAv); %%steepest descent
alpha=linemin(rootparams,presgrad, alpharoot,  MaxPop, Nbins,Tframes,hist, N, Nfac, NexpAv,presgrad); %%finding an adecuate step size
params=rootparams-alpha*presgrad; %%steepest descent update
params(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
params(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
%%

pastconjugdir=-presgrad; %%updating the previous conjugate direction
pastconjugdir(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
pastconjugdir(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
pastparams=rootparams; %%updating the previous argument of the loglikelihood
grad=logligrad(params ,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv); %%gradient of the function to be minimized

%%now we use nonlinear conjugate gradients to fin the MLE
tolerance=10E-5; %%when the gradient's magnitude is smaller than this, the search stops
counter=0;
tic
while sqrt((params-pastparams)'*(params-pastparams))/(MaxPop+1+Nbins)>tolerance %algorithm stops when parameters change on average to the 5th decimal
		   
   conjugdir=getconjugdir(params,pastparams,pastconjugdir, MaxPop,Nbins,Tframes,hist, N, Nfac, NexpAv,grad);  %%preconditioner used is PR, also it can get automatic resets to steepest descent
   conjugdir(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   conjugdir(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
		   
	alpha=linemin(params,conjugdir,alpharoot,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv, grad); %%finding an adecuate step size
		   
   pastparams=params; %%updating the previous argument of the loglikelihood
   pastparams(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   pastparams(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
		   
   params=params+alpha*conjugdir; %%conjugate gradients update
   params(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
   params(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
		   
   pastconjugdir=conjugdir;  %%updating the previous conjugate direction
   %sqrt(logligrad(params ,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv)'*logligrad(params ,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv)) %%checking the norm of the gradient
   grad=logligrad(params ,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv); %%gradient of the function to be minimized
   counter=counter+1;
end
toc
		 
counter

if gauge==0
    Vmle=params(MaxPop+2:end); %%extracting vexation at the minimum
    fmle=params(1:MaxPop+1); %%extracting frustration at the minimum
    [stderrors,CovMatmle]=getCovMat(fmle,Vmle, MaxPop,Nbins,Tframes,0); %getting the covariance matrix for the parameters in the model without gauge transformation 
    fmleError=[0;0;stderrors(1:MaxPop-1)]; %%unpacking errors
    VmleError=stderrors(MaxPop:end); %%unpacking errors
%counter
else
    N=((1:(MaxPop+1))-1)';
    Vmle=params(MaxPop+2:end)-sum(params(MaxPop+2:end))./Nbins; %%extracting vexation at the minimum
    fmle=params(1:MaxPop+1)+sum(params(MaxPop+2:end))*N./Nbins; %%extracting frustration at the minimum
    [stderrors,CovMatmle]=getCovMat(fmle,Vmle, MaxPop,Nbins,Tframes,1);  %getting the covariance matrix for the parameters in the model with gauge transformation that sets the average potential to zero
    fmleError=[0;stderrors(1:MaxPop)]; %%unpacking errors
    VmleError=stderrors(MaxPop+1:end); %%unpacking errors
end



end




