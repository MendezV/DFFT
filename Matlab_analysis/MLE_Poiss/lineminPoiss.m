function alpha=lineminPoiss(params,conjugdir,alpharoot, MaxPop, Nbins, Tframes, hist, N, Nfac, NexpAv,grad)

%%IN
%%-params: Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-conjugdir: a vector of size Nbins  corresponding to the
%%new search direction in the algorithm
%%-alpharoot: scalar that corresponds to the root that will be used in
%%the linesearch algorithm
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


%%performs a backtracking search to determine a sufficiently good step size
%%to guaratee that we are indeed minimizing the function. Sufficiently good
%%is determined by the Armijo-Goldstein condition

%%OUT
%%-alpha: optimal step size along the search direction according to the Armijo-Goldstein condition


tau=0.7; %%control parameter should be between 0 and 1
c1=0.0001; %%control parameter should be between 0 and 1
%c2=0.1;  %%control parameter should be between 0 and 1
m=conjugdir'*grad/sqrt(conjugdir'*conjugdir); %%local slope of the function of alpha  along the search direction
t=-c1*m; %parameter that serves as lower bound in the Armijo-Goldstein condition 
alpha=alpharoot; %%initializing alpha for the linesearch
pastLogli=logliPoiss(params,  MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv);

diff=pastLogli-logliPoiss(params+alpha*conjugdir,MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv); %%difference that will be minimized

%%updating alpha until the Armijo-Goldstein and the sstrong Wolfe conditions are fulfilled, then we
%%are satisfied with the value of alpha that will multiply the conjugate
%%direction in the algorithm

while diff<alpha*t
	alpha=alpha*tau;
    diff=pastLogli-logliPoiss(params+alpha*conjugdir,MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv);
    
end

end

