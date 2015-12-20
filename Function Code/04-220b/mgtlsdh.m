function m = mgtlsdh(d,w,dh)
% MGTLSDH - Global Total Least Squares misfit computation.
% MGTLSDH(D,W,DH) = || sqrtm(W) * (D - DH) ||_F
%
% See GTLS for the possible formats of W.

[sd,N] = size(d);
if length(w(:)) == sd % EWGTLS
  w  = w(:);
  sw = sqrt(w);
  m  = norm(sw(:,ones(N,1)) .* (d - dh),'fro');
else
  sw = chol(w);
  m  = norm(sw * (d - dh),'fro');
end