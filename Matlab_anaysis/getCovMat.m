function [stderrors,CovMat]=getCovMat(f,V, MaxPop,Nbins,Tframes,gauge)

%%IN
%%-f: a MaxPop+1(MaxPop is the maximum number of flies that were observed inside a bin) sized vector that corresponds to the frustration that came out of the MLE including the gauge fixed values
%%-V: a Nbins(number of bins in the system) sized vector that corresponds to the vexation that came out of the MLE
%%-MaxPop:  maximum observed packing in the system
%%-Nbins: total number of bins
%%-Tframes: number of frames
%%-gauge: numerical value, if equal to zero a gauge transformation was
%%performed to the parameters and we have to transform the covariance
%%matrix accordingly


%%Calculates the asymptotic Covariance Matrix as the inverse of the fisher
%%information matrix for our log-likelihood function

%%OUT
%%-stderrors: a (MaxPop-1+Nbins)x1 vector that corresponds to the diagonal of the covariance matrix, corresponding to the
%%variances for each of the parameters that we are estimating without the
%%gauge fixed values, which have no uncertainty
%%-CovMat: a (MaxPop-1+Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators. (if there was no gauge fix)
%%-CovMat:a (MaxPop+Nbins) square, positive, symmetric, invertible
%%Matrix that corresponds to the covariance matrix of the asymptotic
%%gaussian distribution for the ML estimators after performing the gauge transformation
%%which ammounts to performing a similarity transformation to the covriance matrix
%%also we append to the asymptotic covariance the error and covariances of the parameter f(1) that was previously fixed but now has error. (if there was a gauge fix)


%parameters that will be used in the calculation
N=((1:(MaxPop+1))-1)'; %vector with possible occupation numbers in the system

%%frustration sector of the covariance matrix
z=sum(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2); %%size Nbinsx1 normalization for probaility in each bin for our model
probmat=(exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'))./(z*ones(1,MaxPop+1)); %size NbinsxMaxPop+1
HessFF=Tframes*(diag(sum(probmat,1)')-probmat'*probmat); %%this hessian is the fisher information matrix 
%%HessFF=HessFF(3:end,3:end);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case
HessFF=HessFF;

%%vexation sector of the covariance matrix
NensAv=sum((ones(Nbins,1)*N').*exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %% the ensemble average according to our model
NsqensAv=sum((ones(Nbins,1)*((N.^2)')).*exp(-V*N'-ones(Nbins,1)*f')./gamma(ones(Nbins,1)*(N+1)'),2)./z; %% the ensemble average of N^2 according to our model
HessVV=Tframes*diag(NsqensAv-NensAv.*NensAv);

%%mixed sector of the covariance matrix
HessVF=Tframes*probmat.*(ones(Nbins,1)*N'- NensAv*ones(1,MaxPop+1));
%%HessVF=HessVF(:,3:end);%wierd index to eliminate values that dont have uncertainty in the gauge fixed case
HessVF=HessVF;

%%final result
if gauge==0
    
    Hess=[HessFF,HessVF';HessVF,HessVV]; %%hessian, equal to the fisher information matrix
    d=eig(Hess);
    zer=size(diag(d),1)-nnz(single(diag(d).^5)); %% taking out machine precision numbers in the fisher information matrix
    HessFF=HessFF(1+zer:end,1+zer:end);
    HessVF=HessVF(:,1+zer:end);
    Hess=[HessFF,HessVF';HessVF,HessVV]; %%hessian, equal to the fisher information matrix
    if zer==2
        CovMat=inv(Hess); %%pseudoinverse, apparently its singular due to gauge invariance
    else
        extra=zer-2;
        CovMatPre=inv(Hess); %%pseudoinverse, apparently its singular due to gauge invariance
        CovMat=[zeros(extra,extra),zeros(extra,size(CovMatPre,1));zeros(size(CovMatPre,1),extra),CovMatPre];
    end
    stderrors=sqrt(diag(CovMat)); %%the factor of two comes in the taylor expansion for the gaussian approximation
  
else
    Hess=[HessFF,HessVF';HessVF,HessVV]; %%hessian, equal to the fisher information matrix
    d=eig(Hess);
    zer=size(diag(d),1)-nnz(single(diag(d).^5)); %% taking out machine precision numbers in the fisher information matrix
    HessFF=HessFF(1+zer:end,1+zer:end);
    HessVF=HessVF(:,1+zer:end);
    Hess=[HessFF,HessVF';HessVF,HessVV]; %%hessian, equal to the fisher information matrix
    gaugeTMat=eye(size(Hess))+[zeros(size(HessFF)),(N(1+zer:end))*ones(1,Nbins)./Nbins;zeros(size(HessVF)),-ones(size(HessVV))./Nbins];
    CovMatp=(inv(Hess)); %%covariance before gauge transformation as the (well defined) inverse of the  fisher information matrix
    CovMatg=gaugeTMat*CovMatp*gaugeTMat'; %%the gauge transformation to set the average potential is just a linear transformation on the vector of parameters, this is the induced transformation on the covariance matrix
    covf1primevbprime=sum(CovMatp(MaxPop:end,MaxPop:end),2)/(Nbins)+sum(sum(CovMatp(MaxPop:end,MaxPop:end),1),2)/(Nbins^2);
    covf1primefnprime=sum(CovMatp(MaxPop:end,1:MaxPop-1),1)'/(Nbins)+N(3:end)*sum(sum(CovMatp(MaxPop:end,MaxPop:end),1),2)/(Nbins^2);
    redvectcorr=[covf1primefnprime;covf1primevbprime]; %%sector that technically doesnt belong to the asymptotic covariance matrix but its included for having a practical packing of the error propagation on the value of f(1)
    if zer==2
        CovMat=[sum(sum(CovMatp(MaxPop:end,MaxPop:end),1),2)/(Nbins^2),redvectcorr';redvectcorr,CovMatg]; %% covariance matrix with the f(1) uncertainties and covariances added to it, f(1) variance includes a sum of al the variances and coariances of the v's divided by the number of bins squared, the other new sectors are the covariances with the other parameters, which can be easily derived by considering the linearity of the average and the definition of the covariance
    else
        extra=zer-2;
        CovMatPre=[sum(sum(CovMatp(MaxPop:end,MaxPop:end),1),2)/(Nbins^2),redvectcorr';redvectcorr,CovMatg]; %% covariance matrix with the f(1) uncertainties and covariances added to it, f(1) variance includes a sum of al the variances and coariances of the v's divided by the number of bins squared, the other new sectors are the covariances with the other parameters, which can be easily derived by considering the linearity of the average and the definition of the covariance
        CovMat=[zeros(extra,extra),zeros(extra,size(CovMatPre,1));zeros(size(CovMatPre,1),extra),CovMatPre];
    end
    stderrors=sqrt(diag(CovMat)); %diagonal of the matrix above corresponds to the variances of the parameters

end


end
