load Circpred;

counts=occ_f';
Corr;
gauge=0;


MaxPop=max(max(counts)); %maximum observed packing in the system
Nbins=size(counts,1); %total number of bins
Tframes=size(counts,2)/tau; %%number of independent frames (in reality should be downsapled untli tau is equal to 1)

alpharoot =  0.0010000;
rootparams=rand([MaxPop+1+Nbins,1]);

delta=zeros(MaxPop+1,Nbins); %%allocating the matrix of kronecker delta functions
for n=0:MaxPop
delta(n+1,:)=mean(counts'==n); %%basically computing the histogram of counts
				  end
				  hist=delta';
