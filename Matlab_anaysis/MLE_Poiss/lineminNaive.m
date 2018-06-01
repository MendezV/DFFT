function alpha=lineminNaive(params,conjugdir,alpharoot, MaxPop, Nbins, Tframes, hist)

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

%%performs a backtracking search to determine a sufficiently good step size
%%to guaratee that we are indeed minimizing the function. Sufficiently good
%%is determined by the Armijo-Goldstein condition

%%OUT
%%-alpha: optimal step size along the search direction according to the Armijo-Goldstein condition

tau=0.7; %%control parameter should be between 0 and 1
c=0.005; %%control parameter should be between 0 and 1
grad=logligradNaive(params ,  MaxPop, Nbins, Tframes, hist); %%gradient of the function to be minimized
m=conjugdir'*grad/sqrt(conjugdir'*conjugdir); %%local slope of the function of alpha  along the search direction
t=-c*m; %parameter that serves as lower bound in the Armijo?Goldstein condition 
alpha=alpharoot; %%initializing alpha for the linesearch

diff=logliNaive(params,  MaxPop, Nbins, Tframes, hist)-logliNaive(params+alpha*conjugdir,  MaxPop, Nbins, Tframes, hist); %%difference that will be minimized 

%%updating alpha until the Armijo-Goldstein condition is fulfilled, then we
%%are satisfied with the value of alpha that will multiply the conjugate
%%direction in the algorithm
while diff<alpha*t
    alpha=alpha*tau;
    diff=logliNaive(params,  MaxPop, Nbins, Tframes, hist)-logliNaive(params+alpha*conjugdir,  MaxPop, Nbins, Tframes, hist);
end

end

