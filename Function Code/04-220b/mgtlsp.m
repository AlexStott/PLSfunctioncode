function [m,dh] = mgtlsp(d,w,p)
% MGTLSP - Global Total Least Squares misfit computation.
% [M,DH] = MGTLSP(D,W,P) gives the GTLS misfit M and the GTLS
% approximation DH of the data D by the model image(P).
%
% See GTLS for the possible formats of W.

% The computation is based on transformation to a modified TLS   
% problem. As an alternative, see Corollary 2.13.

[sd,N] = size(d);

% Modify P and D
if length(w(:)) == sd % EWGTLS
  w  = w(:); % make it a column vector
  sw = sqrt(w);
  p  = sw(:,ones(1,size(p,2))) .* p;
  d  = sw(:,ones(1,N)) .* d; 
else % GTLS
  sw = chol(w);
  p  = sw * p;
  d  = sw * d;
end

% Compute the TLS misfit for the modified P and D
if nargout == 1
  m = mtlsp(d,p);
else
  [m,dh] = mtlsp(d,p);
end