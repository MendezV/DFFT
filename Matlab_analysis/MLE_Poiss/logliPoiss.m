function logp=logliPoiss(params ,  MaxPop, Nbins, Tframes, hist, N, Nfac, NexpAv)

%%IN
%%-params: Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames
%%-hist: Nbins x (MaxPop+1) matrix in which each row corresponds to the
%%histrogram of counts within each bin
%%-N:  vector of size MaxPop+1 with ordered integers ranging from 0 to MaxPop 
%%-Nfac: vector of size MaxPop+1 with the factorial of ordered integers ranging from 0 to MaxPop 
%%-NexpAv: vector of size Nbins with the average number of individuals within
%%each bin


%%calculates the value of the likelihood function for a given set of
%%parameters and data for our specific model

%%OUT
%%-logp: value of the likelihood function for the parameters and data from
%%the analytic formula in the vexation-only model

%parameters that will be used in the calculation
V=params; %%extracting vexation

%%elements that will be needed for the calculation of the negative log of
%%the probability
z=sum(exp(-V*N')./Nfac,2); %%vector that contains the partition function of each bin, size Nbinsx1
loggammaav=hist*log(gamma(N+1)); %%not necessary for the minimization, could do it but sums to
%a huge number 
logp=Tframes*sum(log(z)+V.*NexpAv+loggammaav)+sum(sum(gammaln(Tframes*hist+1)))-Nbins*gammaln(Tframes+1);


end