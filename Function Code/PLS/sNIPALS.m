function [ Xp, Yp, Xs, Ys, B,W , C,P,b] = sNIPALS( X, Y, n )
%NIPALS PLS algorithm

P  = [];%Matrix for loadings of X onto latent variable T
T = [];%Matrix for latent variables 
Q = [];%Matrix for loadings of Y onto latent variable U
U = [];%Matrix for Y latent variable - taken from T
W = [];%Matrix for weight vector for X
C = [];%Matrix for weight vector for Y
b = [];

Xd = X;
Yd = Y;


for i = 1:n
    u = Y(:,1);
    %for j = 1:10000
    told = 1;
    t = 0;
    while norm(told-t) > 0.001
        w = Xd'*u/(u'*u);
        w = w/norm(w);
        told = t;
        t = Xd*w; %X compute scores
        c = Yd'*t/(t'*t);
        %c = c/norm(c);        
        u = Yd*c;%/(c'*c); %Compute Y scores
        p = Xd'*t/(t'*t); %compute normalised loadings for X
        %p = p/norm(p);
        q = Yd'*u/(u'*u); %compute normalised loadings for Y
        %q = q/norm(q);
        b1 = u'*t/(t'*t);
    end

    
    
    Xd = Xd - t*p';
    Yd = Yd - t*c';
%     Yd = Yd - u*q';
    
    P(:,i) = p;
    Q(:,i) = q;
    T(:,i) = t;
    U(:,i) = u;
    W(:,i) = w;
    C(:,i) = c;
    b(:,i) = b1;
    
    
end

B = (W*inv(P'*W)*C');%W*C';

Xp = P;
Yp = Q;
Xs = T;
Ys = U;

end
