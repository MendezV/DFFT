function grad=logligradPoiss(params , MaxPop, Nbins, Tframes, hist, N, Nfac, NexpAv)

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


%%calculates the gradient of the likelihood function for our specific model
%%from the analytical formula at the value of the parameters being passed as input with a give set of
%%data

%%OUT
%%-grad: a vector of size Nbins that corresponds to the gradient of the
%%likelihood function for the naive model



%parameters that will be used in the calculation
V=params; %%extracting vexation

%Vexation sector of the gradient
z=sum(exp(-V*N')./Nfac,2); %%size Nbinsx1 normalization for probaility in our model
NensAv=sum((ones(Nbins,1)*N').*exp(-V*N')./Nfac,2)./z; %%ensemble average of our model
Vgrad=Tframes*(NexpAv-NensAv); %V sector of the gradient is the difference between the observed average and the model average

grad=Vgrad;

end
