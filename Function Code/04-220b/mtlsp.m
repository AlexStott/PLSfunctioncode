function [m,dh] = mtlsp(d,p)
% MTLSP - Total Least Squares misfit computation.
% [M,DH] = MTLSP(D,P) gives the TLS misfit M and the TLS
% approximation DH of the data D by the model image(P).

% Compute the data approximation dh = p/(p'*p)*p'*d, see Cor 2.13,
% however, because p might be non-minimal, it is safer to use pinv.
dh = p * pinv(p) * d;
m  = norm(d-dh,'fro');