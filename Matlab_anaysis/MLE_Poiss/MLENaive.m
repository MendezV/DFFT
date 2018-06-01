function [VNaive,CovMatNaive,stderrorsNaive]=MLENaive(rootparams,alpharoot,counts,tau)
%%%%%%%%%%
%%IN
%%-rootparams: Nbins values for the vexation at each bin for the naive model, 
%%it corresponds to the root of the minimization algorithm 
%%-rootalpha: scalar that corresponds to the root that will be used in
%%the linesearch algorithm that is implemented in the minimization 
%%-counts: a NbinsxTframes matrix with the number of flies observed in each
%%bin at each timeframe
%%-tau:correlation time corrects for the number of independent time samples


%%finds the minimum of g(V)=-logP(V|data) by the method of preconditioned conjugate
%%gradients the argument of this function at the minimum gives the maximum
%%likelihood estimates for the parameters in the naive model. 

%%OUT
%%-VNaive: a Nbins(number of bins in the system) sized vector that corresponds
%%to the vexation that came out of the MLE 
%%-stderrors: a (MaxPop-1+Nbins)x1 vector that corresponds to the diagonal of the covariance matrix, corresponding to the
%%variances for each of the parameters that we are estimating without the
%%gauge fixed values, which have no uncertainty
%%-CovMatNaive: a (Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators for the naive model. 





%%%%%%%%
MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)


delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
    delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
end
hist=delta';

%%first iteration is steepest descent with a variable step size 

presgrad=-logligradNaive(rootparams , MaxPop, Nbins, Tframes, hist); %%steepest descent
alpha=lineminNaive(rootparams,presgrad,alpharoot, MaxPop, Nbins, Tframes, hist); %%finding an adecuate step size
params=rootparams+alpha*presgrad; %%steepest descent update
%%

pastconjugdir=presgrad; %%updating the previous conjugate direction
pastparams=rootparams; %%updating the previous argument of the loglikelihood

%%now we use nonlinear conjugate gradients to fin the MLE
tolerance=0.01; %%when the gradient's magnitude is smaller than this, the search stops
%counter=0;
while sqrt(logligradNaive(params , MaxPop, Nbins, Tframes, hist)'*logligradNaive(params ,  MaxPop, Nbins, Tframes, hist))>tolerance
    
   conjugdir=getconjugdirNaive(params,pastparams,pastconjugdir, MaxPop, Nbins, Tframes, hist);  %%preconditioner used is PR, also it can get automatic resets to steepest descent
   alpha=lineminNaive(params,conjugdir,alpharoot, MaxPop, Nbins, Tframes, hist) %%finding an adecuate step size
   pastparams=params;  %%updating the previous argument of the loglikelihood
   params=params+alpha*conjugdir; %%conjugate gradients update
   pastconjugdir=conjugdir;  %%updating the previous conjugate direction
   sqrt(logligradNaive(params , MaxPop, Nbins, Tframes, hist)'*logligradNaive(params , MaxPop, Nbins, Tframes, hist)) %%checking the norm of the gradient
   %counter=counter+1
end

VNaive=params; %%extracting vexation at the minimum
[stderrorsNaive,CovMatNaive]=getCovMatNaive(VNaive, MaxPop, Nbins, Tframes); %getting the covariance matrix for the parameters in the model without gauge transformation 
%counter

end

