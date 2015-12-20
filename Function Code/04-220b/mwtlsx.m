function [m,dh] = mwtlsx(d,w,x)
% MWTLSX - Weighted Total Least Squares misfit computation.
% [M,DH] = MWTLSX(D,W,X) gives the WTLS misfit M and the WTLS
% approximation DH of the data D by the model B(X).
%
% See the help of WTLS for the possible formats of W.

[m,p] = size(x);
if nargout == 1
  if m > p
    m = mwtlsr(d,w,x2r(x));
  else
    m = mwtlsp(d,w,x2p(x));
  end
else  
  if m > p
    [m,dh] = mwtlsr(d,w,x2r(x));
  else
    [m,dh] = mwtlsp(d,w,x2p(x));
  end
end
