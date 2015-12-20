function r = minr(r)
% MINR - Minimal kernel representation.
% MINR(R) is a minimal kernel representation of ker(R).
% MINR(R,TOL), TOL is a user defined tolerance for 
% computing numerical rank.

if nargin == 1
  tol = 1e-14; % default 
end
[p,sd] = size(r);

[u,s,v] = svd(r,'econ');
s       = diag(s); % column vector
pmin    = sum(s > tol);
if pmin < p % truncate R to a minimal one
  r = s(1:pmin,ones(1,sd)) .* v(:,1:pmin)';
end
