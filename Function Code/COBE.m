function [ Abar ] = COBE( Y, n, eps )

Q = zeros(size(Y,1),size(Y,1),n);
H = zeros(size(Y,1),size(Y,2),n);

for i = 1:n
    Yn = Y(:,:,i);
    [U, E, V] = svd(Yn);
    Q(:,:,i) = U;
    H(:,:,i) = E*V';
end

f = 0;

% z = zeros(size(Y,1),n);
z = Y(:,1:n);
zprev = 10*ones(size(Y,1),n);
aprev = 10*ones(size(Y,1),1);
abar = zeros(size(Y,1),1);
Abar = [];
k = 0;

while f<eps
    aprev = 10*ones(size(Y,1),1);
    abar = zeros(size(Y,1),1);
    iternum = 0;
    while norm(abar - aprev) > 0.001
        aprev = abar;
        abar = zeros(size(Y,1),1);
        for i = 1:n
            Qmult = Q(:,:,i);            
            abar = abar + Qmult*z(:,i);
        end
        abar = abar/norm(abar,'fro');
        zprev = z;
        for i = 1:n
            Qmult = Q(:,:,i);
            z(:,i) = Qmult*abar;
        end
        iternum = iternum + 1;
        if iternum > 10000
            return 
        end
    end
    iternum
    f = 0;
    for i = 1:n
        Qmult = Q(:,:,i);
        f=f+norm(Qmult*z(:,i)-abar,'fro')^2;
    end
    f
    Abar = [Abar abar];
    for i = 1:n
        Q(:,:,i) = Q(:,:,i)*(eye(size(Y,1)) - z(:,i)*z(:,i)');
    end
    if k > size(Y,2)
        return
    end
    k = k+1
end

            

end

