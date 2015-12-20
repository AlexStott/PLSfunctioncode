function [M,dh] = mwtlsp(d,w,p)
% MWTLSP - Weighted Total Least Squares misfit computation.
% [M,DH] = MWTLSP(D,W,P) gives the WTLS misfit M and the WTLS
% approximation DH of the data D by the model image(P).
%
% See the help of WTLS for the possible formats of W.

% See Theorem 2.12
[sd,N] = size(d);
dh = zeros(sd,N);
switch size(w,2)
 case N % EWTLS
  m  = size(p,2);    % # of inputs
  for i = 1:N
    tmp = w(:,i*ones(1,m)) .* p;
    dh(:,i) = p / ( p' * tmp ) * tmp' * d(:,i);
  end
  M = mwtlsdh(d,w,dh);
 case sd % WTLS
  M  = 0;
  for i = 1:N
    tmp     = w(:,:,i) * p;
    dh(:,i) = p / ( p' * tmp ) * tmp' * d(:,i);
    ddi     = d(:,i) - dh(:,i);
    M       = M + ddi' * w(:,:,i) * ddi;
  end
  M = sqrt(M);
 case sd*N % FWTLS
  p   = kron(eye(N),p);
  tmp = w * p;
  dh  = p / ( p' * tmp ) * tmp' * d(:);
  dd  = d(:) - dh;
  M   = sqrt(dd' * w * dd);
  dh  = reshape(dh,sd,N);
end