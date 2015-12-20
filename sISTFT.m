function [ x ] = sISTFT( X,w,over,N )

fftlen = length(w);

numffts = size(X,2);

x = zeros(N+fftlen,1);

for i = 0:numffts-1
    Temp = ifft([X(:,i+1); conj(X(end-1:-1:2,i+1))])./w.';
    x(1+(i)*(fftlen-over):(i+1)*(fftlen-over)) = real(Temp(1:(fftlen-over)));
end

x = x(1:N-fftlen);


end

