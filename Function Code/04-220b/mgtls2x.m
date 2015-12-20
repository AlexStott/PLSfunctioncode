function m = mgtlsx(d,wl,wr,x)
% MGTLS2X - Global Total Least Squares misfit computation.
% [M,DH] = MGTLS2X(D,WL,WR,X) gives the GTLS2 misfit M and the 
% GTLS2 approximation DH of the data D by the model B(X).
%
% See GTLS for the possible formats of WL and WR.

[m,p] = size(x);
if nargout == 1
  if m > p
    m = mgtls2r(d,wl,wr,x2r(x));
  else
    m = mgtls2p(d,wl,wr,x2p(x));
  end
else  
  if m > p
    [m,dh] = mgtls2r(d,wl,wr,x2r(x));
  else
    [m,dh] = mgtls2p(d,wl,wr,x2p(x));
  end
end