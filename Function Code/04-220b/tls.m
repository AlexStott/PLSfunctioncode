function [r,p,M,dh] = tls(d,m)
% TLS - Total Least Squares approximation.
%
% [r,p,M,dh] = tls(d,m)
%
% D = [d1 ... dN] - data matrix
% m  - complexity specification, m < size(D,1)
% R  - parameter of a kernel representation of the TLS model 
% P  - parameter of an image representation of the TLS model
% M  - TLS misfit 
% DH - TLS data approximation 

% The algorithm is an application of Theorem 2.8.

sd = size(d,1);

% If DH is not needed, more efficient is first 
% to compress the data by QR, see Note 9.
if nargout < 4 
  d = triu(qr(d'))';
  d = d(:,1:sd);
end

[u,s,v] = svd(d,'econ');
s = diag(s)'; % row vector
r = u(:,m+1:sd)';
p = u(:,1:m);
M = norm(s(m+1:sd));

if nargout >= 4
  dh = ( u(:,1:m) .* s(ones(sd,1),1:m) ) * v(:,1:m)';
end