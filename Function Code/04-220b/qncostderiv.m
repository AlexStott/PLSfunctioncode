% QNCOSTDERIV - WTLS cost function and first derivative
% evaluation for the function FMINUNC.

function [M,dM] = qncostderiv(x,d,v,m)

% constants
[sd,N] = size(d);
p      = sd - m;
m1     = m + 1;
deriv  = (nargout == 2);
if size(v,3) == 1
  c = 1;
else
  c = 2;
end

e  = x' * d(1:m,:) - d(m1:sd,:);  % residual
M  = 0;
dM = zeros(m,p);
for i = 1:N
  ei  = e(:,i);
  if c == 1
    vax = v(1:m,i*ones(1,p)) .* x;
    yi  = ( x' * vax + diag(v(m1:sd,i)) ) \ ei;
  else
    vax = v(1:m,1:m,i) * x - v(1:m,m1:sd,i);
    vbx = v(m1:sd,1:m,i) * x - v(m1:sd,m1:sd,i);
    yi  = ( x' * vax - vbx ) \ ei;
  end
  M   = M + ei' * yi;
  if deriv
    dM = dM + ( d(1:m,i) - vax * yi ) * yi';
  end  
end
M = sqrt(M);