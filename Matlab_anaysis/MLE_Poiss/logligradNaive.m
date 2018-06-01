function grad=logligradNaive(params , MaxPop, Nbins, Tframes, hist)

%%IN
%%-params: Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames

%%calculates the gradient of the likelihood function for our specific model
%%from the analytical formula at the value of the parameters being passed as input with a give set of
%%data

%%OUT
%%-grad: a vector of size Nbins that corresponds to the gradient of the
%%likelihood function for the naive model



%parameters that will be used in the calculation
V=params; %%extracting vexation
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system


%Vexation sector of the gradient
z=sum(exp(-V*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in our model
NensAv=sum((ones(Nbins,1)*N').*exp(-V*N')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average of our model
NexpAv=hist*N;%observed experimental average number of flies at each bin
Vgrad=Tframes*(NexpAv-NensAv); %V sector of the gradient is the difference between the observed average and the model average

grad=Vgrad;

end
