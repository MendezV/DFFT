function logp=logliNaive(params ,  MaxPop, Nbins, Tframes, hist)

%%IN
%%-params: Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames

%%calculates the value of the likelihood function for a given set of
%%parameters and data for our specific model

%%OUT
%%-logp: value of the likelihood function for the parameters and data from
%%the analytic formula in the naive model

%parameters that will be used in the calculation
V=params; %%extracting vexation
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system

%%elements that will be needed for the calculation of the negative log of
%%the probability
z=sum(exp(-V*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%vector that contains the partition function of each bin, size Nbinsx1
NexpAv=hist*N; %experimental average number of flies at each bin
loggammaav=hist*log(gamma(N+1)); %%not necessary for the minimization, could do it but sums to
%a huge number 
logp=Tframes*sum(log(z)+V.*NexpAv+loggammaav)+sum(sum(gammaln(Tframes*hist+1)))-Nbins*gammaln(Tframes+1);


end