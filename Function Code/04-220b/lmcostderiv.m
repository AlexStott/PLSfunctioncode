% LMCOSTDERIV - WTLS cost function and first derivative
% evaluation for the function LSQNONLIN.

function [M,dM] = lmcostderiv(x,d,v,m)

% constants
[sd,N] = size(d);
p      = sd - m;
m_1    = m + 1;
mp     = m * p;
jac    = (nargout == 2);
wtls   = (size(v,3) > 1);

a    = d(1:m,:)';
b    = d(m_1:sd,:)';
xext = [x;-eye(p)];
e = (a * x - b)';  % residual
M = zeros(p,N);
if jac
  dM = zeros(N*p,m*p);
end
for i = 1:N
  if wtls
    tmp = chol(xext' * v(:,:,i) * xext);
  else
    tmp = chol(xext' * (v(:,i*ones(1,p)) .* xext));
  end
  M(:,i) = tmp' \ e(:,i);
  if jac
    invGi     = inv(xext' * v(:,:,i) * xext);
    invsqrtGi = chol(invGi);
    A     = kron(eye(p),a);
    vaix  = v(1:m,1:m,i) * x;
    yi    = invGi * e(:,i);
    for j = 1:mp
      E   = zeros(m,p); E(j) = 1;
      f   = A((i-1)*p+1:i*p,j);
      dMj = invsqrtGi * f;
      DG  = E' * vaix + vaix' * E - E' * v(1:m,m_1:sd,i) ...
            - v(m_1:sd,1:m,i) * E;
      dMj = dMj - 1/2 * invsqrtGi * DG * yi;
      dM((i-1)*p+1:i*p,j) = dM((i-1)*p+1:i*p,j) + dMj;
    end
  end
end
M = M(:);
