function tau=Corr(counts)

fnmin=0;

% ===============================================================
% Autocorrelation function (using FFT), averaged over boxes
c=(mean( abs(ifft(abs(fft(counts)).^2))/size(counts,1) - ones(size(counts,1),1)*mean(counts).^2 ,2));
%figure(1); plot(c(1:100),'rx')

lit=[1:5]'; % first few points
p=polyfit(lit,log(c(lit)),1); 
fprintf('Correlation time (in frames)...\n');
tau=-1/p(1)


end