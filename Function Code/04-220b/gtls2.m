function [r,p,M,dh] = gtls2(d,wl,wr,m,tol)
% GTLS2 - Global Total Least Squares approximation 
% with two side weighting.
%
% [r,p,M,dh] = gtls2(d,wl,wr,m,tol)
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% WL - positive definite sd x sd (left) weight matrix or an sd x 1
%      vector wl, such that WL = diag(wl) (element-wise weighting)
% WR - positive definite N x N (right) weight matrix or an N x 1
% m  - complexity specification, m < sd 
% TOL - tolerance for checking ill conditioning (default 1e-14)
%      vector wr, such that WR = diag(wr) (element-wise weighting)
% R  - parameter of a kernel representation of the GTLS2 model 
% P  - parameter of an image representation of the GTLS2 model
% M  - GTLS2 misfit 
% DH - GTLS2 data approximation 

% The algorithm is an application of Theorem 2.9.

[sd,N] = size(d); 

% GTLS2 or EW-GTLS2 case?
if length(wl(:)) == sd
  wl = wl(:);  % make it a column vector
  wr = wr(:)'; % make it a row vector
  c = 1;
else
  c = 0;
end

% Check the conditioning of WL and WR
if nargin == 4
  tol = 1e-14;
end
if c
  if any(wl < tol)
    error('Ill conditioned weight matrix WL.')
  end
  if any(wr < tol)
    error('Ill conditioned weight matrix WR.')
  end
else
  if rcond(wl) < tol
    error('Ill conditioned weight matrix WL.')
  end
  if rcond(wr) < tol
    error('Ill conditioned weight matrix WR.')
  end
end

% Modified data
if c
  swl = sqrt(wl);
  swr = sqrt(wr);
  d  = swl(:,ones(1,N)) .* d .* swr(ones(sd,1),:);  
else
  swl = chol(wl);
  swr = chol(wr);
  d   = swl * d * swr;
end

% Find the TLS approximation
if nargout < 4
  [r,p,M] = tls(d,m);
else
  [r,p,M,dh] = tls(d,m);
end

% Transform back
if c
  r = r .* swl(:,ones(1,sd-m))';
  if nargout >= 2
    p = swl(:,ones(1,m)) .\ p;
    if nargout >= 4
      dh = swl(:,ones(1,N)) .\ dh ./ swr(ones(sd,1),:);
    end
  end  
else
  r = r * swl;
  if nargout >= 2
    p = swl \ p;
    if nargout >= 4
      dh = swl \ dh / swr;
    end
  end
end