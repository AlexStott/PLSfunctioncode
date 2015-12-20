function [x,M,dh] = wtls(d,w,m)
% WTLS - Weighted Total Least Squares approximation
%
% [x,M,dh] = wtls(d,w,m)
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% W  - sd x N matrix (EWTLS case)
%      sd x sd x N tensor, such that Wi = W(:,:,i) (WTLS case),
%      sd.N x sd.N weight matrix for vec(D) (FWTLS case)
% m  - complexity specifification, m < sd
% X  - parameter of an I/O representation of the WTLS model
% M  - WTLS misfit
% DH - WTLS data approximation 
%
% Note: requires Optimization Toolbox.

% Determine the case
[sd,N] = size(d);
switch size(w,2)
 case N % EWTLS
  c = 1;
 case sd % WTLS
  c = 2;
 case sd*N % FWTLS
  error('FWTLS case not implemented yet.');
 otherwise
  error('Wrong dimension of W.')
end

x0  = r2x(wtlsini(d,w,m)); % initial approximation
w   = w2v(w,c);            % done once for all iterations
opt = optimset('GradObj','on','LargeScale','off','Display','off');
[x,M,flag,out] = fminunc(@(x)costderiv(x, d, w), x0, opt);
if nargout > 2
  [M,dh] = mwtlsx(d,w,x);
else
  M = sqrt(M);
end

% ----------------------------------

function [M,dM] = costderiv(x,d,v)
% COSTDERIV - WTLS cost function and first derivative evaluation.

% constants
[sd,N] = size(d);
[m,p]  = size(x);
m1     = m + 1;
deriv  = (nargout == 2);
wtls   = (size(v,3) > 1);

e  = x' * d(1:m,:) - d(m1:sd,:);  % residual
M  = 0;
dM = zeros(m,p);
% Recognize is it WTLS or EWTLS problem
for i = 1:N
  ei  = e(:,i);
  if wtls
    vax = v(1:m,1:m,i) * x - v(1:m,m1:sd,i);
    vbx = v(m1:sd,1:m,i) * x - v(m1:sd,m1:sd,i);
    yi  = ( x' * vax - vbx ) \ ei;
  else
    vax = v(1:m,i*ones(1,p)) .* x;
    yi  = ( x' * vax + diag(v(m1:sd,i)) ) \ ei;
  end
  M   = M + ei' * yi;
  if deriv
    dM = dM + ( d(1:m,i) - vax * yi ) * yi';
  end  
end

