function p = minp(p,tol)
% MINP - Minimal image representation.
% MINP(P) is a minimal image representation of image(P).
% MINP(P,TOL), TOL is a user defined tolerance for 
% computing numerical rank.

if nargin == 1
  tol = 1e-14; % default 
end
[sd,m] = size(p);

[u,s] = svd(p,'econ');
s       = diag(s)'; % row vector
mmin    = sum(s > tol);
if mmin < m % truncate P to a minimal one
  p = u(:,1:mmin) .* s(ones(sd,1),1:mmin);
end