function [p,info,dh] = wtlsap(d,w,m,opt)
% WTLSAP - Weighted Total Least Squares approximation
% by alternating projections (Algorithm 2.1).
%
% [p,info,dh] = wtlsap(d,w,m,opt)
%
% D = [d1 ... dN] - data matrix, sd := size(D,1)
% W - sd x N weight matrix (EWTLS problem), 
%     sd x sd x N tensor (WTLS problem), or 
%     sd.N x sd.N weight matrix for vec(d) (FWTLS)
% m    - complexity specification, m < sd
% OPT  - options for the optimization algorithm, see OPTIMSET
%   OPT.MAXITER - maximum number of iterations 
%   OPT.TOLFUN  - convergence tolerance for the function value  
%   OPT.P0      - user defined initial approximation 
% P    - parameter of an image representation of the WTLS model
% INFO - structure containing exit information:
%   INFO.M    - WTLS misfit
%   INFO.TIME - execution time
%   INFO.ITER - number of iterations performed
%   Note: INFO.ITER = OPT.MAXITER indicates lack of convergence
% DH   - WTLS data approximation 

tic % measure the execution time

% Constants
[sd,N] = size(d);
sdm    = sd*m;

% Which case?
switch size(w,2)
 case N
  c = 1;
 case sd
  c = 2;
 case sd*N
  c = 3;
end

% Default parameters
if nargin > 3
  if isempty(opt.MaxIter)
    opt.MaxIter = 1000;
  end
  if isempty(opt.TolFun)
    opt.TolFun = 1e-10;
  end
  % Compute initial approximation if not given
  if isfield(opt,'p0') 
    p = opt.p0;
    [M,dh] = mwtlsp(d,w,p);  
  elseif isfield(opt,'P0') 
    p = opt.P0;
    [M,dh] = mwtlsp(d,w,p);
  else
    [r,p,M,dh] = wtlsini(d,w,m);
  end
else
  opt = optimset('Display','off');  
  opt.TolFun  = 1e-10;
  opt.MaxIter = 1000;
  [r,p,M,dh]  = wtlsini(d,w,m);
end
l = zeros(m,N); % reserve memory for L

% Main iteration loop
cont = 1;
iter = 0;
while (cont)
  % Solve the relaxation problems with a special method
  % depending on the weight matrix structure 
  if c == 1 % EWTLS
    % Solve RLX1
    for i = 1:N
      wip = w(:,i*ones(1,m)) .* p;
      l(:,i) = (p'*wip)\wip'*d(:,i);
    end
    % Solve RLX2
    A  = zeros(sdm);
    b  = zeros(sdm,1);
    for k = 1:N
      wk  = w(:,k);
      for i = 1:m
        for j = 1:m
          A((i-1)*sd+1:i*sd,(j-1)*sd+1:j*sd) = diag( ...
              diag(A((i-1)*sd+1:i*sd,(j-1)*sd+1:j*sd)) ...
              + l(i,k) * l(j,k) * wk);
        end
        b((i-1)*sd+1:i*sd) = b((i-1)*sd+1:i*sd) + ...
            l(i,k) * (wk .* d(:,k));
      end
    end
    p  = reshape(A\b,sd,m);    
   elseif c == 2 % WTLS
    % Solve RLX1
    for i = 1:N
      ptwi = p' * w(:,:,i); 
      l(:,i) = (ptwi*p)\ptwi*d(:,i);
    end
    % Solve RLX2
    A = zeros(sdm);
    b = zeros(sdm,1);
    for k = 1:N
      wk = w(:,:,k);
      for i = 1:m
        for j = 1:m
          A((i-1)*sd+1:i*sd,(j-1)*sd+1:j*sd) = ...
              A((i-1)*sd+1:i*sd,(j-1)*sd+1:j*sd) + ...
              l(i,k) * l(j,k) * wk;
        end
        b((i-1)*sd+1:i*sd) = b((i-1)*sd+1:i*sd) + ...
            l(i,k) * wk * d(:,k);
      end
    end
    p  = reshape(A\b,sd,m);    
   else % FWTLS
    % Solve RLX1
    bp = kron(eye(N),p);
    vl = (bp'*w*bp)\bp'*w*d(:);
    l  = reshape(vl,m,N);
    % Solve RLX2
    bl = kron(l',eye(sd));
    vp = (bl'*w*bl)\bl'*w*d(:);
    p  = reshape(vp,sd,m);
  end
  
  % Display
  iter  = iter + 1;
  Mold  = M;
  dh    = p*l;
  M     = mwtlsdh(d,w,dh);
  switch lower(opt.Display)
   case 'iter'
    fprintf('%2d : Misfit = %18.16f\n',k,M)
  end  
  
  % Exit condition
  rerr  = abs(M-Mold)/M;
  cont = (iter < opt.MaxIter) & (rerr > opt.TolFun); 
end

info.M    = M;
info.iter = iter;
info.time = toc;
