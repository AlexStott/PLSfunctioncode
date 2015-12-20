function [x,info,dh] = wtlsopt(d,w,m,opt)
% WTLSOPT - Weighted Total Least Squares approximation
% by standard local optimization algorithms (Algorithm 2.4).
%
% [x,info,dh] = wtlsopt(d,w,m,opt)
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% W - sd x N weight matrix (EWTLS problem), 
%     sd x sd x N tensor (WTLS problem), or 
%     sd.N x sd.N weight matrix for vec(d) (FWTLS)
% m - complexity specification, m < sd
% OPT - options for the optimization algorithm, see OPTIMSET
%   OPT.X0     - user defined initial approximation
%   OPT.ALG    - optimization algorithm 
%                (fminunc, lsqnonlin, or fminsearch)
% X    - parameter of an I/O representation of the WTLS model
% INFO - structure containing exit information:
%   INFO.M     - WTLS misfit
%   INFO.TIME  - execution time
%   INFO.ITER  - number of iterations performed
%   INFO.FEVAL - number of cost function evaluations
%   Note: INFO.ITER = OPT.MAXITER indicates lack of convergence
% DH   - WTLS data approximation 
%
% Note: requires Optimization Toolbox.

tic % measure the execution time

% Default optimization method
if nargin > 3
  if isfield(opt,'alg')
    alg = opt.alg;
  elseif isfield(opt,'ALG') 
    alg = opt.ALG;
  elseif isfield(opt,'Alg') 
    alg = opt.Alg;
  else % default
    alg = 'fminunc'; 
  end
  if isempty(opt.Display)
    opt = optimset(opt,'Display','off');
  end
  % Default initial approximation
  if isfield(opt,'x0') 
    x0 = opt.x0;
  elseif isfield(opt,'X0') 
    x0 = opt.X0;
  else % default
    x0 = r2x(wtlsini(d,w,m));
  end
else % default
  opt = optimset('Display','off');
  alg = 'fminunc'; 
  x0 = r2x(wtlsini(d,w,m));
end

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

% Run the optimization method
switch lower(alg)
 case 'fminunc'
  opt = optimset(opt,'GradObj','on','LargeScale','off');
  v = w2v(w,c); % done once for all iterations
  [x,M,flag,out] = fminunc(@(x)qncostderiv(x,d,v,m),x0,opt);
  M = sqrt(M);
 case 'lsqnonlin'
  opt = optimset(opt,'Jacobian','off','LevenbergMarquardt', ...
                 'on','LargeScale','off','DerivativeCheck','off');
  v = w2v(w); % done once for all iterations
  [x,M,res,flag,out] = lsqnonlin(@(x)lmcostderiv(x,d,v,m), ...
                                 x0, [], [], opt);
  M = sqrt(M);
 case 'fminsearch'  
  [x,M,flag,out] = fminsearch(@(x)mwtlsx(d,w,x), x0, opt);
 otherwise
  error([alg ' optimization method not supported.']);
end
  
% Assign output variables
info.M     = M;
info.iter  = out.iterations;
info.feval = out.funcCount;
info.time  = toc;
if nargout > 2
  [M,dh] = mwtlsx(d,w,x);
end
