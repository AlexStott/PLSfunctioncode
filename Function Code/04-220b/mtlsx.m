function [m,dh] = mtlsx(d,x)
% MTLSX - Total Least Squares misfit computation.
% [M,DH] = MTLSX(D,X) gives the TLS misfit M and the TLS
% approximation DH of the data D by the model B(X).

if nargout == 1
  if m > p
    m = mtlsr(d,x2r(x));
  else
    m = mtlsp(d,x2p(x));
  end
else  
  if m > p
    [m,dh] = mtlsr(d,x2r(x));
  else
    [m,dh] = mtlsp(d,x2p(x));
  end
end
