function [ X ] = sSTFT( x, w, over )

N = length(x);
fftlen = length(w);

xapp = [x zeros(1,fftlen)];
numffts = ceil(N/(fftlen-over));
X = zeros(1+fftlen/2,numffts);

for i = 1:numffts
    Temp = fft(w.*xapp(1+(i-1)*(fftlen-over):(i-1)*(fftlen-over)+fftlen));
    X(:,i) = Temp(1:1+fftlen/2).';
end
    
X = X(:,1:numffts-1);

end

