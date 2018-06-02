function [VPoiss,CovMatPoiss,VPoissError]=MLEPoiss(root,alpharoot,counts,tau)
%%%%%%%%%%
%%IN
%%-root: Nbins values for the vexation at each bin for the naive model, 
%%it corresponds to the root of the minimization algorithm, if equal to
%%'random' the initial guess for the parameter vector has entries uniformly
%%distributed between 0 and one
%%-rootalpha: scalar that corresponds to the root that will be used in
%%the linesearch algorithm that is implemented in the minimization 
%%-counts: a NbinsxTframes matrix with the number of flies observed in each
%%bin at each timeframe
%%-tau:correlation time corrects for the number of independent time samples


%%finds the minimum of g(V)=-logP(V|data) by the method of preconditioned conjugate
%%gradients the argument of this function at the minimum gives the maximum
%%likelihood estimates for the parameters in the Poisson model. 

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


%%
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
%%
%setting up initial guess for the parameter vector
if root=='random'
   rootparamsP=rand([Nbins,1]);
end
%%first iteration is steepest descent with a variable step size 
%%
presgrad=logligradPoiss(rootparamsP , MaxPop, Nbins, Tframes, hist, N, Nfac, NexpAv); %%steepest descent
alpha=lineminPoiss(rootparamsP,-presgrad,alpharoot, MaxPop, Nbins, Tframes, hist, N, Nfac, NexpAv,presgrad); %%finding an adecuate step size
params=rootparamsP-alpha*presgrad; %%steepest descent update


pastconjugdir=-presgrad; %%updating the previous conjugate direction
pastparams=rootparamsP; %%updating the previous argument of the loglikelihood
grad=logligradPoiss(params ,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv); %%gradient of the function to be minimized
%%
%%now we use nonlinear conjugate gradients to fin the MLE
tolerance=10E-7; %%when the gradient's magnitude is smaller than this, the search stops
counter=0;

tic
while sqrt((params-pastparams)'*(params-pastparams))/(MaxPop+1+Nbins)>tolerance %algorithm stops when parameters change on average to the 5th decimal
		   
   conjugdir=getconjugdirPoiss(params,pastparams,pastconjugdir, MaxPop,Nbins,Tframes,hist, N, Nfac, NexpAv,grad);  %%preconditioner used is PR, also it can get automatic resets to steepest descent	   
   alpha=lineminPoiss(params,conjugdir,alpharoot,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv, grad); %%finding an adecuate step size	   
   pastparams=params; %%updating the previous argument of the loglikelihood
   params=params+alpha*conjugdir; %%conjugate gradients update   
   pastconjugdir=conjugdir;  %%updating the previous conjugate direction
   grad=logligradPoiss(params , MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv); %%gradient of the function to be minimized
   counter=counter+1;
end
toc
		 
counter
%%

VPoiss=params; %%extracting vexation at the minimum
[VPoissError,CovMatPoiss]=getCovMatPoiss(VPoiss, MaxPop, Nbins, Tframes); %getting the covariance matrix for the parameters in the model without gauge transformation 
%counter

end

