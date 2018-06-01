function conjugdir=getconjugdirNaive(params,pastparams,pastconjugdir, MaxPop, Nbins, Tframes, hist)

%%IN
%%-params: Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-pastparams corresponds to the position for past
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-pastconjugdir: a vector of size Nbins  corresponding to the
%%last search direction in the algorithm
%%-counts: a NbinsxTframes matrix with the number of flies observed in each
%%bin at each timeframe


%%Gets the search direction in the non-linear conjugate gradients minimization algorithm by implementing a preconditioner given by the
%%Polak-Ribiere (Fletcher-Reeves) formula with a modification to generate automatic resets to
%%steepest descent if necesary
%%


%%OUT
%%-conjugdir: a vector of size Nbins  corresponding to the
%%new search direction in the algorithm

pastgrad=logligradNaive(pastparams , MaxPop, Nbins, Tframes, hist); %gradient in the last iteration
grad=logligradNaive(params , MaxPop, Nbins, Tframes, hist); %gradient in this iteration
betaPR=grad'*(grad-pastgrad)/(pastgrad'*pastgrad); %Polak-Ribiere 
beta=max([0,betaPR]); %%automatic reset for steepest decent
%betaFR=grad'*grad/(pastgrad'*pastgrad); %Fletcher-Reeves
%beta=betaFR;
conjugdir=-grad+beta*pastconjugdir; %%conjugate direction in the algorithm
end