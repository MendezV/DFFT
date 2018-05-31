
n=counts';
fnmin=0;

% ===============================================================
% Autocorrelation function (using FFT), averaged over boxes
c=(mean( abs(ifft(abs(fft(n)).^2))/size(n,1) - ones(size(n,1),1)*mean(n).^2 ,2));
figure(1); plot(c(1:100),'rx')

lit=[1:5]'; % Yunus data seems to show two time scales??
p=polyfit(lit,log(c(lit)),1); 
fprintf("Correlation time (in frames)...\n");
tau=-1/p(1)
fprintf("\n");
