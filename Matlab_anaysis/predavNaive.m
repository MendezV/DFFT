function [predictedAvNaive, errorbarsNaive, muNaive]=predavNaive(VNaive,CovMatNaive,Nflies)

%%IN:
%%-VNaive: a Nbins(number of bins in the system) sized vector that
%%corresponds to the vexation that came out of the MLE for the naive model
%%-CovMatNaive a (MaxPop-1+Nbins) square covariance matrix obtained from
%%getCovMat
%%-Nflies: Number of flies in the sistem that we are trying to predict
%%-MaxPop: the maximu number of flies that can occupsiy a bin in the system
%%we are trying to predict

%%predicts the average number of flies in each bin using the Naive model and
%%then determines the uncertainties taking into account correlations
%%between the parameters in the model

%%OUT:
%%-predictedAv a vector of size Nbins (the number of bins in the system) in
%%which each entry is the average number of flies as predicted from the
%%Naive
%%model implemented with a binary search that determines the chemical
%%potential mu that fixes the total number of flies (vector is sorted in the same order as the occupations matrix first dimension)
%%-errorbars a vector of size Nbins (the number of bins in the system) in
%%which each entry is the one sigma error on the average number of flies as
%%predicted from the Naive
%%model for each bin(vector is sorted in the same order as the occupations matrix first dimension)
%%-muNaive:Chemical otential that fixes the average total umber of flies in
%%the system, for the Naive model without interactions

MaxPop=100;
Nbins=size(VNaive,1); %total number of bins
N=((1:(MaxPop+1))-1)';


%% first we make the prediction 
NmaxSteps=130;    %number of steps before giving up on the search
tolerance= 1e-7;  %tolerance in the difference between the actual number of flies and the one that comes up with the guessed mu
counter=0;
mu=max(VNaive); %%convenient root for the search algorithm
NfliesGuess=0.0;
dist=abs(max(VNaive)-min(VNaive));
a=mu+2*dist;  %guess for the highest lower bound for the chemical potential
b=mu-2*dist;   %guess for the lowest higher bound for the chemical potential
remu=0;

while counter<NmaxSteps
   
    mu=(a+b)/2;
    z=sum(exp(-(VNaive-mu)*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in our model
    predictedAvNaive=sum((ones(Nbins,1)*N').*exp(-(VNaive-mu)*N')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average in the naive model
    NfliesGuess=sum(predictedAvNaive); %% total number of flies for this iteration of mu
    
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
z=sum(exp(-(VNaive-remu)*N')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in each bin for our model
NsqensAv=sum((ones(Nbins,1)*(N.^2)').*exp(-(VNaive-remu)*N')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average of n^2 in our model
NensAv=sum((ones(Nbins,1)*N').*exp(-(VNaive-remu)*N')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %%ensemble average of our model
partialV=-(NsqensAv-NensAv.*NensAv);
vexerror=sqrt(diag(CovMatNaive)); %%errors taken from the diagonal of the matrix CovMatVNaive

%%correlations between values of V do not appear as the partial derivative with respect to
%%vexations in other bins vanishes, the average in each bin doesnt depend on the
%%vexation in other bins

errorbarsNaive=sqrt((partialV.^2).*(vexerror.^2));

muNaive=remu;  %%output the chemical potential
end