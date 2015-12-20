function [ P, Q, T, U, BETA  ] = ssimpls( X, Y, ncomp )

S = X.'*Y;
T  = zeros(size(X,1),ncomp);
U  = zeros(size(Y,1),ncomp);
P = zeros(size(X,2),ncomp);
Q = zeros(size(Y,2),ncomp);
V = P;
R = [];

for i = 1:ncomp
    
    qmat = eig(S);
    qvec = qmat(:,1);
    r = S*qvec;
    t = X*r;
    r = r/sqrt(t.'*t);
    t = t/sqrt(t.'*t);
    p = X.'*t;
    q = Y.'*t;
    u = Y*q;
    v=p;
    if i > 1
        v = v - V*(V.'*p);
        u = u - T*(T.'*u);
    end
    v = v/sqrt(v.'*v);
    S = S - v*(v.'*S);
    
    T(:,i) = t;
    U(:,i) = u;
    Q(:,i) = q;
    P(:,i) = p;
    V(:,i) = v;
    R = [R r];

end

BETA = R*Q.';
    
end

