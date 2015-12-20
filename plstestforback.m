clear all
close all
clc
N = 1000;
f1 = 0.01;
t = [1:N];
f2 = 0.06;
x1 = randn(N,1);%-sin(f1.*t)';%+1i*sin(f2.*t)';
x2 = 5*randn(N,1) + 0.2*x1;%sign(x1);

x3 =  0.2*x1+ 5*randn(N,1);%sin(f2.*t)';

x4 = 0.01*randn(N,1);
x5 = 0.1*randn(N,1);

Xe = [x1 x2 x3];% x2+0.3*x3];% x4-mean(x4)];
Xa = [zeros(N,1) zeros(N,1) x4];
Ya = [zeros(N,1) 0.2*x5 x5];
% Mix = [1 2+1i 0+0.03i;2 0.2-0.2i 2; 1+1i 3 0.2-1i];
% Mix = [1 2 0 1;4 2 0.2 2;-1 1 3 0.2;1 2 3 -0.2];
Mix = [2 0.3 2;1 0.2 3;1 1 3];%[4; 2; 3];%
% Mix = [2 4 ; 0  1;1 1 ];
Ye = Xe*Mix;
Xbef = Xe;
% X = X + [x4 x4 x4];
% Y = Y + [x5 0.1*x5 x5];
Noise = [0.1*randn(N,3)];
X = Xe + Noise;
% X(:,3) = X(:,3)+x5;
Y =Ye +  0.1*randn(N,3);
% Y(:,1) = Y(:,1) + 1*randn(N,1);
% Y(:,2) = Y(:,2) + sin(f2.*t)';

[Pf, Qf, Tf, Uf, Bnf, Wf, Cf] = snplssvdunitsc(X,Y,3);
[Pb, Qb, Tb, Ub, Bnb, Wb, Cb] = snplssvdunitsc(Y,X,3);

Xg = Tf*(Pf.'*Wf)*Wf.';
Xginv = Wf*inv(Pf.'*Wf)*Tf.';
Xmpinv = pinv(Xg);

Yg = Tb*(Pb.'*Wb)*Wb.';
Yginv = Wb*inv(Pb.'*Wb)*Tb.';

Bnf*Bnb*Bnf
Bnf
Bnb*Bnf*Bnb
Bnb
Bnf*Bnb
Bnb*Bnf

norm(Bnf*Bnb*Bnf-Bnf,'fro')
norm(Bnb*Bnf*Bnb-Bnb,'fro')

% [Pf2, Qf2, Tf2, Uf2, Bnf2, Wf2, Cf2] = snplssvd(X,Y,3);
% 
% [Ux, Ex, Vx] = svd(X);
% [Uy, Ey, Vy] = svd(Y);
% Tx = Ux*Ex;
% Tx = Tx(:,1:3);
% Px = Vx(:,1:3);
% 
% Xtilpc = Tx*Px';
% 
% Ty = Uy*Ey;
% Ty = Ty(:,1:3);
% Py = Vy(:,1:3);
% 
% Ytilpc = Ty*Py';
% 
% Cxpcr = Tx\Y;
% 
% Yrecpcr = Tx*Cxpcr;
% 
% Bpcr = (Px)*Cxpcr;


% Cypcr = Ty\X;

% [Pb2, Qb2, Tb2, Ub2, Bnb2, Wb2, Cb2] = snplssvd(Y,X,3);

% Xg2 = Tf2*(Pf2.'*Wf2)*Wf2.';
% Xginv2 = Wf2*inv(Pf2.'*Wf2)*Tf2.';

% [Pf, Qf, Tf, Uf, Bf, Rf, Vf] = ssimpls(X,Y,3);
% [Pb, Qb, Tb, Ub, Bb, Rb, Vb] = ssimpls(Y,X,3);



% Rfopt = Rf;
% Rbopt = Rb;
% 
% for i = 1:10
% Rfopt = Qb*Rb'*Qf*pinv(Qf'*Rb*(Qb'*Qb)*Rb'*Qf);
% h = Qb*Rbopt'*Qf;
% H =h'*h;
% % gradRf = 2*Rfopt*H-2*h;
% 
% % Rfopt = Rfopt - 0.01*gradRf;
% 
% Tfopt = X*Rfopt;
% Rfopt(:,1) = Rfopt(:,1)/norm(Tfopt(:,1));
% Rfopt(:,2) = Rfopt(:,2)/norm(Tfopt(:,2));
% % Rfopt(:,3) = Rfopt(:,3)/norm(Tfopt(:,3));
% Tfopt(:,1) = Tfopt(:,1)/norm(Tfopt(:,1));
% Tfopt(:,2) = Tfopt(:,2)/norm(Tfopt(:,2));
% % Tfopt(:,3) = Tfopt(:,3)/norm(Tfopt(:,3));
% Pfopt = X'*Tfopt;
% Qfopt = Y'*Tfopt;
% Ufopt = Y*Qfopt;
% Bfopt = Rfopt*Qfopt';
% 
% Qf = Qfopt;
% 
% % Rbopt = pinv(Bfopt'*Bfopt)*(Bfopt'*Qb)*pinv(Qb'*Qb);
% 
% % gradRb = 2*Bfopt'*Bfopt*Rbopt*Qb'*Qb - 2*Bfopt'*Qb;
% 
% % Rbopt = Rbopt - 0.01*gradRb;
% 
% % Tbopt = Y*Rbopt;
% % Rbopt(:,1) = Rbopt(:,1)/norm(Tbopt(:,1));
% % Rbopt(:,2) = Rbopt(:,2)/norm(Tbopt(:,2));
% % % Rbopt(:,3) = Rbopt(:,3)/norm(Tbopt(:,3));
% % Tbopt(:,1) = Tbopt(:,1)/norm(Tbopt(:,1));
% % Tbopt(:,2) = Tbopt(:,2)/norm(Tbopt(:,2));
% % % Tbopt(:,3) = Tbopt(:,3)/norm(Tbopt(:,3));
% % Pbopt = Y'*Tbopt;
% % Qbopt = X'*Tbopt;
% % Ubopt = X*Qbopt;
% Bbopt = Rbopt*Qb';
% 
% % Qb=Qbopt;
% % Rb = Rbopt;
% 
% norm(Bfopt*Bbopt-eye(3),'fro')
% 
% end
% 
% 
% norm(Bfopt*Bbopt-eye(3),'fro')
% norm(Bf*Bb-eye(3),'fro')
% 
% norm((Y - X*Bfopt),'fro')
% norm((Y - X*Bf),'fro')
% norm((X - Y*Bbopt),'fro')
% norm((X - Y*Bb),'fro')
