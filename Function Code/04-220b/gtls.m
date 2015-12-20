function [r,p,M,dh] = gtls(d,w,m,tol)
% GTLS - Global Total Least Squares approximation with 
% one side weighting.
%
% [r,p,M,dh] = gtls(d,w,m,tol)
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% W  - positive definite sd x sd weight matrix or an sd x 1
%      vector w, such that W = diag(w) (element-wise weighting)
% m  - complexity specification, m < sd
% TOL - tolerance for checking ill conditioning (default 1e-14)
% R  - parameter of a kernel representation of the GTLS model 
% P  - parameter of an image representation of the GTLS model
% M  - GTLS misfit 
% DH - GTLS data approximation 

% The algorithm is an application of Theorem 2.9

[sd,N] = size(d); 

% GTLS or EW-GTLS case?
if length(w(:)) == sd
  w = w(:); % make it a column vector
  c = 1;
else
  c = 0;
end

% Check the conditioning of W
if nargin == 3
  tol = 1e-14; % default
end
if c
  if any(w < tol)
    error('Ill conditioned weight matrix W.')
  end
else
  if rcond(w) < tol
    error('Ill conditioned weight matrix W.')
  end
end

% Modified data
if c
  sw = sqrt(w);
  d  = sw(:,ones(1,N)) .* d;
else
  sw = chol(w);
  d  = sw * d;
end

% Find the TLS approximation
if nargout < 4
  [r,p,M] = tls(d,m);
else
  [r,p,M,dh] = tls(d,m);
end

% Transform back
if c
  r = r .* sw(:,ones(1,sd-m))';
  if nargout >= 2
    p = sw(:,ones(1,m)) .\ p;
    if nargout >= 4  
      dh = sw(:,ones(1,N)) .\ dh;
    end
  end 
else
  r = r * sw;
  if nargout >= 2
    p = sw \ p;
    if nargout >= 4  
      dh = sw \ dh;
    end
  end
end
