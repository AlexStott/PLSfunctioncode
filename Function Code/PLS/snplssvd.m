function [ Xp, Yp, Xs, Ys, B,W , C,b] = snplssvd( X, Y, n )
%NIPALS PLS algorithm computed from SVD http://www.eigenvector.com/evriblog/?p=86

P  = [];%Matrix for loadings of X onto latent variable T
T = [];%Matrix for latent variables 
Q = [];%Matrix for loadings of Y onto latent variable U
U = [];%Matrix for Y latent variable - taken from T
W = [];%Matrix for weight vector for X
C = [];%Matrix for weight vector for Y
b= [];

S = X'*Y;

for i = 1:n
    [Us,Es,Vs] = svd(S); %obtain eigen value for S
    w1 = Us(:,1); %Dominant Eigen vector is the weight 
    t1 = X*w1; %X compute scores
    %t1 = t1/norm(t1);
    c1 = Y'*t1/(t1'*t1); %compute Y loading for latent variable
    u1 = Y*c1; %Compute Y scores
   % u1 = u1/norm(u1);
    p1 = X'*t1/(t1'*t1); %compute normalised loadings for X
    q1 = Y'*u1/(u1'*u1); %compute normalised loadings for Y
    b1 = u1'*t1/(t1'*t1);
    
    X = X - t1*p1';
    Y = Y - t1*c1';
    
    P(:,i) = p1;
    Q(:,i) = q1;
    T(:,i) = t1;
    U(:,i) = u1;
    W(:,i) = w1;
    C(:,i) = c1;
    b(:,i) = b1;
    
    S = X'*Y;
end

B = (W*inv(P'*W)*C');

Xp = P;
Yp = Q;
Xs = T;
Ys = U;

end
