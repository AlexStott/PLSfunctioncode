function [x,info,dh] = wtlspr(d,w,m,opt)
% WTLSPR - Weighted Total Least Squares approximation
% by the algorithm of Premoli-Rastelo (Algorithm 2.2).
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% W    - sd x N matrix (EWTLS problem), 
%     sd x sd x N tensor (WTLS problem), or 
%     sd.N x sd.N weight matrix for vec(d) (FWTLS case)
% m    - complexity specification, m < sd
% OPT  - options for the optimization algorithm, see OPTIMSET
%   OPT.X0    - user defined initial approximation 
% X    - parameter of an I/O representation of the WTLS model
% INFO - structure containing exit information:
%   INFO.M    - WTLS misfit
%   INFO.TIME - execution time
%   INFO.ITER - number of iterations performed
%   Note: INFO.ITER = OPT.MAXITER indicates lack of convergence
% DH   - WTLS data approximation 

tic % measure the execution time

[sd,N] = size(d);

% Determine the case
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

% Default parameters
if nargin > 3
  if isempty(opt.MaxIter)
    opt.MaxIter = 20;
  end
  if isempty(opt.TolX)
    opt.TolX = 1e-5;
  end
  if isempty(opt.TolFun)
    opt.TolFun = 1e-5;
  end
  if isempty(opt.Display)
    opt.Display = 'off';
  end
  if isfield(opt,'x0') 
    x = opt.x0;
  elseif isfield(opt,'X0') 
    x = opt.X0;
  else
    x = r2x(wtlsini(d,w,m));
  end
else
  opt.MaxIter = 20;
  opt.TolX    = 1e-5;
  opt.TolFun  = 1e-5;
  opt.Display = 'off';
  x = r2x(wtlsini(d,w,m));
end

% Main iteration loop
a    = d(1:m,:)';
b    = d(m+1:sd,:)';
p    = sd - m; % # of outputs
mp   = m * p ; % often used constant 
cont = 1;
k    = 0;
M    = 0;
while (cont)
  e    = (a * x - b)';   % residual matrix
  xext = [x;-eye(p)];    % extended parameter
  g    = zeros(mp);
  h    = zeros(mp,1);
  Mold = M;
  M    = 0;              % misfit
  for i = 1:N
    if c == 1
      mi  = inv(xext' * (w(:,i*ones(1,p)) .* xext));
    else
      mi  = inv(xext' / w(:,:,i) * xext);
    end
    yi = mi * e(:,i);
    ai = a(i,:)';
    if c == 1
      g   = g + kron(mi,ai*ai') - kron(yi*yi',diag(w(1:m,i)));
      tmp = ai * b(i,:) * mi;      
    else
      g   = g + kron(mi,ai*ai') - kron(yi*yi',w(1:m,1:m,i));
      tmp = ai * b(i,:) * mi - w(1:m,m+1:sd,i) * yi*yi';
    end
    h   = h + tmp(:);
    M   = M + e(:,i)' * yi;
  end
  M     = sqrt(M);
  xold  = x;
  x     = reshape(g\h,m,p);
  k     = k + 1;
  dX    = norm(x - xold,'fro') / norm(x,'fro');
  dM    = abs(M - Mold) / M;
  % Display
  switch lower(opt.Display)
   case 'iter'
    fprintf('%2d : Misfit = %18.16f\n',k,M)
  end
  % Exit condition
  cont = (k < opt.MaxIter) & (dX > opt.TolX) & (dM > opt.TolFun);
end

info.M = M;
info.iter = k;
info.time = toc;
if nargout > 2
  [M,dh] = mwtlsx(d,w,x);
end
