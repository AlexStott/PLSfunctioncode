function m = mwtlsdh(d,w,dh)
% MWTLSDH - Weighted Total Least Squares misfit computation.
% [M,DH] = MWTLSDH(D,W,DH) gives the WTLS misfit M achived by
% the data approximation DH.
%
% See the help of WTLS for the possible formats of W.

[sd,N] = size(d);
dd = d - dh;
switch size(w,2)
 case N % EWTLS
  m = norm( sqrt(w) .* dd, 'fro');
 case sd % WTLS
  m  = 0;
  for i = 1:N
    m = m + dd(:,i)' * w(:,:,i) * dd(:,i);
  end
  m = sqrt(m);  
 case sd*N % FWTLS
  m = sqrt( dd(:)' * w * dd(:) );
end