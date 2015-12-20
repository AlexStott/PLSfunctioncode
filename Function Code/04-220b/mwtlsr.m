function [M,dh] = mwtlsr(d,w,r)
% MWTLSR - Weighted Total Least Squares misfit computation.
% [M,DH] = MWTLSR(D,W,R) gives the WTLS misfit M and the WTLS
% approximation DH of the data D by the model ker(R).
%
% See the help of WTLS for the possible formats of W.

% See Theorem 2.10
[sd,N] = size(d);
dd = zeros(sd,N);
switch size(w,2)
 case N % EWTLS
  p = size(r,1); % # of putputs
  for i = 1:N
    tmp     = w(:,i*ones(1,p)) .\ r';
    dd(:,i) = tmp / ( r * tmp ) * r * d(:,i);
  end
  M = norm( sqrt(w) .* dd, 'fro');
 case sd % WTLS
  M  = 0;
  for i = 1:N
    tmp     = w(:,:,i) \ r';
    dd(:,i) = tmp / ( r * tmp ) * r * d(:,i);
    M       = M + dd(:,i)' * w(:,:,i) * dd(:,i);
  end
  M = sqrt(M);
 case sd*N % FWTLS
  r   = kron(eye(N),r);
  tmp = w \ r';
  dd  = tmp / ( r * tmp ) * r * d(:);
  M   = sqrt(dd' * w * dd);
  if nargout > 1
    dd = reshape(dd,sd,N);
  end
end  
if nargout > 1
  dh = d - dd;
end
