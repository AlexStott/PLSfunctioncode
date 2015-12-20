function [ x ] = invcwt( X,nv )

dt = 1/1200; s0 = 2*dt; ds = 0.4875; 
% NbSc = 32;
as = s0*2.^((0:nv-1)*ds)';

sca = repmat(as,[1,size(X,2)]);
sqsca = sqrt(sca);

Wr = real(X); 
x = sum(Wr./sqrt(sqsca),1);


N = size(X,2);
for i = 1:length(as);
    w=as(i)*linspace(0,N-1,N)*2*pi/N;
    mu=2*pi;
    cmu = (1+exp(-mu^2)-2*exp(-3/4*mu^2)).^(-1/2);
    kmu = exp(-1/2*mu^2);
    morlet =cmu*pi^(-1/4)*(exp(-1/2*(mu-w).^2)-kmu*exp(-1/2*w.^2));
%     morlet(i,:)=sqrt(as(i))*morlet;
    morlet(i,:)=sqrt(as(i))*morlet;
end



Wdelta = sum(morlet,2)/N;
RealWdelta = real(Wdelta);
RealWdelta = RealWdelta(:);
C = sum(RealWdelta./sqrt(sca(:,1)));
% x = (1/C)*x;

tab_MUL = [...
    9.3414	7.7459	6.3586	5.1910	4.2406	3.4910	2.9152 ...
    2.4802	2.1528	1.9062	1.7165	1.5668	1.4467	1.3499 ...
    1.2730	1.2116	1.1611	1.1163	1.0762	1.0362	1.0014 ...
    0.9820	0.9799	0.9971	1.0332	1.0748	1.1006	1.0956 ...
    1.0550	1.0016	0.8997	0.8060	0.7218	0.6961	0.7392 ...
    0.8171	0.9192	1.0350	1.1008	1.1464 ...
    ];
tab_VAL = (1:0.25:10);
param = 6;
D = tab_VAL-param;
[mini,idx] = min(abs(D));

if mini<sqrt(eps)
    mulWAV = tab_MUL(idx);
else
    I1 = find(D<0,1,'last');
    I2 = find(D>0,1,'first');
    T1 = tab_VAL(I1);
    T2 = tab_VAL(I2);
    mulWAV = ((T2-param)*tab_MUL(I1)+  ...
        (param-T1)*tab_MUL(I2))/(T2-T1);
end

% x = (x-mean(x))/mulWAV;



end
