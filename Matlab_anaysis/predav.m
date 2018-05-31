function [predictedAv, errorbars, muDFT]=predav(f,V,CovMatF,CovMatV,Nflies,sameexperiment,gauge)

%%IN:
%%-f: a MaxPop+1(MaxPop is the maximum number of flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values
%%-V: a Nbins(number of bins in the system) sized vector that corresponds to the vexation that came out of the MLE
%%-CovMatV a (MaxPop-1+Nbins) square covariance matrix obtained from
%%getCovMatF: it has less dimensions than expected from the number of parameters in the model due to the lack of error
%%in the gauge fixed parameters.
%%-CovMatF a (MaxPop-1+Nbins) square covariance matrix obtained from
%%getCovMatF: it has less dimensions than expected from the number of parameters in the model due to the lack of error
%%in the gauge fixed parameters.
%%-Nflies: Number of flies in the sistem that we are trying to predict
%%-sameexperiment: int, if it is 1 we take into account correlations between parameters V and F otherwise they will not be correlated 
%%-gauge numerical value, if 0 no gauge transformation is applied to the
%%parameters and f(0)=f(1)=0. else the gauge is set so that the average
%%potential is zero and f(1) corresponds to the sum of the previous values
%%for the potential


%%predicts the average number of flies in each bin using the DFT model and
%%then determines the uncertainties taking into account correlations
%%between the parameters in the model

%%OUT:
%%-predictedAv a vector of size Nbins (the number of bins in the system) in
%%which each entry is the average number of flies as predicted from the DFT
%%model implemented with a binary search that determines the chemical
%%potential mu that fixes the total number of flies (vector is sorted in the same order as the occupations matrix first dimension)
%%-errorbars a vector of size Nbins (the number of bins in the system) in
%%which each entry is the one sigma error on the average number of flies as predicted from the DFT
%%model for each bin(vector is sorted in the same order as the occupations matrix first dimension)
%%-muDFT:Chemical otential that fixes the average total umber of flies in
%%the system, for the DFT model
MaxPop=size(f,1)-1; %maximum observed packing in the system where F was extracted (must be greater or equal to the maximum observed packing in the data we are about to predict)
Nbins=size(V,1); %total number of bins in the data we are about to predict
N=((1:(MaxPop+1))-1)'; %%possible bin occupations for the dataset that we are trying to predict


%% first we make the prediction 

NmaxSteps=130;    %number of steps before giving up on the search
tolerance= 1e-7;  %tolerance in the difference between the actual number of flies and the one that comes up with the guessed mu
counter=0;
mu=max(V); %%convenient root for the search algorithm
NfliesGuess=0.0;
dist=abs(max(V)-min(V));
a=mu+2*dist;  %guess for the highest lower bound for the chemical potential
b=mu-2*dist;   %guess for the lowest higher bound for the chemical potential
remu=0;

while counter<NmaxSteps
   
    mu=(a+b)/2;
    z=sum(exp(-(V-mu)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in our model
    predictedAv=sum((ones(Nbins,1)*N').*exp(-(V-mu)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average of our model
    NfliesGuess=sum(predictedAv); %% total number of flies for this iteration of mu
    
   if abs(NfliesGuess-Nflies)<tolerance %% if the number of flies is right we have the "real mu" and break the while
       remu=mu; %%real mu in case we need it or something
       break
   end
   
   counter=counter+1; %% if the number of flies is not right counter goes up and we modify accordingly for the binary search
   if(NfliesGuess<Nflies)
		b=mu;
		
   else
		a=mu;
   end
    
end


%% then we calculate the error on that prediction by standard propagation of errors
%%these procedurres were implemented in a similar way to calculate elements
%%of the covariance matrix but we have to repeat them due to the presence
%%of the chemical potential we have just introduced

%%with respect to V
z=sum(exp(-(V-remu)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in each bin for our model
NsqensAv=sum((ones(Nbins,1)*(N.^2)').*exp(-(V-remu)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average of n^2 in our model
NensAv=sum((ones(Nbins,1)*N').*exp(-(V-remu)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average of our model
partialV=-(NsqensAv-NensAv.*NensAv);
stderrors=sqrt(diag(CovMatV)); %%errors taken from the diagonal of the matrix CovMatV due to the separability of the distriution the covariances between vexations don't enter in the calculation
vexerror=stderrors(end-Nbins+1:end); %%unpacking this way we dont need to determine the gauge


%%correlations between values of V do not appear as the partial derivative with respect to
%%vexations in other bins vanishes, the average in each bin doesnt depend on the
%%vexation in other bins

%%with respect to f
probmat=(exp(-(V-remu)*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'))./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1
partialF=-probmat.*(ones(Nbins,1)*N'- NensAv*ones(1,MaxPop+1));

if gauge==0 %condition if the gauge is fixed so that the frustration is set to be zero for f(0) and f(1)
   
    partialF=partialF(:,3:end); %wierd index to eliminate values that dont have uncertainty in the gauge fixed case
    
    %%correlations if F do appear in the calculation as the ensemble average in each bin depends
    %%in all the f parameters

    partialFpartialF=zeros(Nbins,MaxPop-1,MaxPop-1);
    corrFF=CovMatF(1:MaxPop-1,1:MaxPop-1);
    ferrorsq=zeros(Nbins,1);

    for i=1:Nbins
        partialFpartialF(i,:,:)=partialF(i,:)*partialF(i,:)';
        ferrorsq(i)=sum(sum(squeeze(partialFpartialF(i,:,:)).*corrFF,2),1);
    end

else %%condition if the gauge is fixed so that the average potential is zero
    
    partialF=partialF(:,2:end); %wierd index to eliminate values that dont have uncertainty in the gauge fixed case

    %%correlations if F do appear in the calculation as the ensemble average in each bin depends
    %%in all the f parameters, but we also have to take into account the
    %%correlations with the value of f(1)

    partialFpartialF=zeros(Nbins,MaxPop,MaxPop);
    corrFF=CovMatF(1:MaxPop,1:MaxPop); %%multiplying off diagonal elemnts by 2 right off the bat
    ferrorsq=zeros(Nbins,1);

    for i=1:Nbins
        partialFpartialF(i,:,:)=partialF(i,:)*partialF(i,:)';
        ferrorsq(i)=sum(sum(squeeze(partialFpartialF(i,:,:)).*corrFF,2),1);
    end

end


if gauge==0   %condition if the gauge is fixed so that the frustration is set to be zero for f(0) and f(1)
   
    %%for both taking into account correlated errors
    partialFpartialV=partialF.*(partialV*ones(1,MaxPop-1)); %%formula is the multiplication of the partial derivatives times the covariance
    corrFV=CovMatF(MaxPop:end,1:MaxPop-1); %%add two to each index for the non gauge fixed case
else  %%condition if the gauge is fixed so that the average potential is zero
   
    %%for both taking into account correlated errors and including
    %%correlations with f(1)
    partialFpartialV=partialF.*(partialV*ones(1,MaxPop)); %%formula is the multiplication of the partial derivatives times the covariance
    corrFV=CovMatF(MaxPop+1:end,1:MaxPop); %%add two to each index for the non gauge fixed case

end
 

%%now we calculate the error bars with standard propagation of errors to
%%first order taking into account correlations 

if sameexperiment==1 %%condition if both V and F come from the same experiment
    errorbars=sqrt((partialV.^2).*(vexerror.^2)+ferrorsq+2*sum(partialFpartialV.*corrFV,2));
else  %%condition if both V and F come from different experiments
    errorbars=sqrt((partialV.^2).*(vexerror.^2)+ferrorsq);
end

muDFT=remu; %%output the chemical potential

end