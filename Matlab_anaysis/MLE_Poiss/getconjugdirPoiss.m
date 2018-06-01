function conjugdir=getconjugdirPoiss(params,pastparams,pastconjugdir, MaxPop, Nbins, Tframes, hist, N, Nfac, NexpAv,grad)

%%IN
%%-params: Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-pastparams corresponds to the position for past
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-pastconjugdir: a vector of size Nbins  corresponding to the
%%last search direction in the algorithm
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames
%%-hist: Nbins x (MaxPop+1) matrix in which each row corresponds to the
%%histrogram of counts within each bin
%%-N:  vector of size MaxPop+1 with ordered integers ranging from 0 to MaxPop 
%%-Nfac: vector of size MaxPop+1 with the factorial of ordered integers ranging from 0 to MaxPop 
%%-NexpAv: vector of size Nbins with the average number of individuals within
%%each bin
%%-grad:  vector of size MaxPop+1 + Nbins that corresponds to the gradient of the likelihood function



%%Gets the search direction in the non-linear conjugate gradients minimization algorithm by implementing a preconditioner given by the
%%Polak-Ribiere (Fletcher-Reeves) formula with a modification to generate automatic resets to
%%steepest descent if necesary
%%


%%OUT
%%-conjugdir: a vector of size Nbins  corresponding to the
%%new search direction in the algorithm


pastgrad=logligradPoiss(pastparams ,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv); %gradient in the last iteration

%%Polak-Ribiere 
%betaPR=grad'*(grad-pastgrad)/(pastgrad'*pastgrad); %Polak-Ribiere 
%beta=max([0,betaPR]); %%automatic reset for steepest decent

%%%Fletcher-Reeves
betaFR=grad'*grad/(pastgrad'*pastgrad); %Fletcher-Reeves
beta=betaFR;

%%Hestenes-Stiefel
%betaHS=-grad'*(grad-pastgrad)/(pastconjugdir'*(grad-pastgrad)); %Hestenes-Stiefel
%beta=betaHS;

%%Dai?Yuan
%betaDY=-grad'*grad/(pastconjugdir'*(grad-pastgrad)); %Dai?Yuan
%beta=betaDY;

conjugdir=-grad+beta*pastconjugdir; %%conjugate direction in the algorithm

end