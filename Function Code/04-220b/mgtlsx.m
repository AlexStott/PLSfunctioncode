function [m,dh] = mgtlsx(d,w,x)
% MGTLSX - Global Total Least Squares misfit computation.
% [M,DH] = MGTLSX(D,W,X) gives the GTLS misfit M and the GTLS
% approximation DH of the data D by the model B(X).
%
% See GTLS for the possible formats of W.

[m,p] = size(x);
if nargout == 1
  if m > p
    m = mgtlsr(d,w,x2r(x));
  else
    m = mgtlsp(d,w,x2p(x));
  end
else  
  if m > p
    [m,dh] = mgtlsr(d,w,x2r(x));
  else
    [m,dh] = mgtlsp(d,w,x2p(x));
  end
end