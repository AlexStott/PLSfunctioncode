function [m,dh] = mgtlsr(d,w,r)
% MGTLSR - Global Total Least Squares misfit computation.
% [M,DH] = MGTLSR(D,W,R) gives the GTLS misfit M and the GTLS
% approximation DH of the data D by the model ker(R).
%
% See GTLS for the possible formats of W.

% The computation is based on a transformation to a modified 
% TLS problem. As an alternative, see Corollary 2.11.

[sd,N] = size(d);

% Modify R and D
if length(w(:)) == sd % EWGTLS
  w  = w(:); % make it a column vector
  sw = sqrt(w);
  r  = r ./ sw(:,ones(1,size(r,1)))';
  d  = sw(:,ones(1,N)) .* d; 
else % GTLS
  sw = chol(w);
  r  = r / sw;
  d  = sw * d;
end

% Compute the TLS misfit for the modified R and D
if nargout == 1
  m = mtlsr(d,r);
else
  [m,dh] = mtlsr(d,r);
end