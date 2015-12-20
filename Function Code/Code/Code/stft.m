function  [S]=stft(s,N)
        
%         t = -2:0.001:2;    
% 
% % Start @ 100Hz, cross 200Hz at t=1 sec 
% y = chirp(t,100,1,200,'quadratic')+(randn(1,length(t))*0); 
% s=y;

%%  STFT

%N=500;
M=length(s);
n=M;
t1=fliplr(s);
x=[zeros(1,N/2),s,zeros(1,N/2)]; %zero pad each side in order to preserve the length of the original signal of applying STFT


% hmin=s(1:N);
% hmax=t1(1:N);
% x=[hmin,s,hmax]; %reflecting boundary conditions




n=length(x);

t=1:1:n-N;
S1=zeros(N/2,length(t));
tau=(0:N-1);

w=gausswin(N,20)';    %Guassian window with a variable window size   

 
 for k=1:length(t)
    temp=fft(x(t(k)+tau).*w);
    g(:,k)=ifft(temp)./w;
    S1(:,k)=temp(1:N/2);
end

%S=S1(:,(N/2)+1:(M+(N/2)));

%S=S1(:,(N)+1:(M+(N)));
S=S1;

