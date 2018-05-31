function grad=logligrad(params ,MaxPop,Nbins,Tframes, hist, N, Nfac, NexpAv)

%%IN
%%-params: a vector of size MaxPop+1 + Nbins (MaxPop+1 values for the frustration,the number of flies 
%%can go from zero to the maximum observed packing in all the bins. And
%%Nbins values for the vexation at each bin ) corresponds to the position for current
%%iteration of the non-linear conjugate gradients minimization algorithm
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames


%%calculates the gradient of the likelihood function for our specific model
%%from the analytical formula at the value of the parameters being passed as input with a give set of
%%data

%%OUT
%%-grad: a vector of size MaxPop+1 + Nbins that corresponds to the gradient of the likelihood function

%parameters that will be used in the calculation
V=params(MaxPop+2:end); %%extracting vexation
f=params(1:MaxPop+1); %%extracting frustration
%f(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
%f(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix


%Vexation sector of the gradient
z=sum(exp(-V*N'-ones(Nbins,1)*f')./Nfac,2); %%size Nbinsx1 normalization for probaility in our model
NensAv=sum((ones(Nbins,1)*N').*exp(-V*N'-ones(Nbins,1)*f')./Nfac,2)./z; %%ensemble average of our model
Vgrad=Tframes*(NexpAv-NensAv); %V sector of the gradient is the difference between the observed average and the model average

%frustration sector of the gradient
probmat=(exp(-V*N'-ones(Nbins,1)*f')./Nfac)./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1
fgrad=Tframes*(sum(hist-probmat,1))'; % f sector of the gradient as the difference of the model probability and the hitogram
fgrad(1)=0; %%fixing the gauge to avoid singularity in the covariance matrix
fgrad(2)=0; %%fixing the gauge to avoid singularity in the covariance matrix
%complete gradient
grad=[fgrad;Vgrad];

end
