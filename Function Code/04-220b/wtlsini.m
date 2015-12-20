function [r,p,M,dh] = wtlsini(d,w,m)
% WTLSINI - Initial approximation for the WTLS problem.
%
% [r,p,M,dh] = wtlsini(d,w,m)
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% W - sd x N weight matrix (EWTLS problem), 
%     sd x sd x N tensor (WTLS problem), or 
%     sd.N x sd.N weight matrix for vec(d) (FWTLS)
% m   - complexity specification, m < sd
% R/P - parameters of kernel/image representation of the model
% M   - WTLS misfit achieved by the model 
% DH  - data approximation

% For initial approximation, we use the EWGTLS2 solution with
% WL and WR, such that WL*WR' is the best rank-1 approx.
% of W in the EWTLS case, or the matrix formed from the diagonal 
% elements of W1,...,WN in the WTLS case, or the matrix 
% formed from the diagonal of W in the FWTLS case. See Note 11.

[sd,N] = size(d);

% Check the given weight matrix and choose method
switch size(w,2)
 case N  % EWTLS
  [u,s,v] = svd(w);
 case sd % WTLS
  % extract the diagonals of the weight matrics
  ew = zeros(sd,N);
  for i = 1:N
    ew(:,i) = diag(w(:,:,i));
  end  
  [u,s,v] = svd(ew);
 case sd*N % General  
  ew = reshape(diag(w),sd,N);
  [u,s,v] = svd(ew);  
 otherwise
  error('Wrong dimension of W.')
end

% Best rank-1 approximation of W or WE
s  = sqrt(s(1));
wl = abs(u(:,1) * s);
wr = abs(v(:,1) * s);

% Compute the EWGTLS2 approximation with these weights
[r,p]  = gtls2(d,wl,wr,m);
[M,dh] = mwtlsr(d,w,r);