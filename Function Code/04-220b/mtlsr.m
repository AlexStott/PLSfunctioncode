function [m,dh] = mtlsr(d,r)
% MTLSR - Total Least Squares misfit computation.
% [M,DH] = MTLSR(D,R) gives the TLS misfit M and the TLS
% approximation DH of the data D by the model ker(R).

% Compute the data correction dd = r'*inv(r*r')*r*d, see Cor 2.11,
% however, because r might be non-minimal, it is safer to use pinv.
dd = pinv(r) * r * d; 
m  = norm(dd,'fro');
if nargout > 1
  dh = d - dd;
end