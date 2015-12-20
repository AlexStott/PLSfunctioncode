function m = mgtls2dh(d,wl,wr,dh)
% MGTLS2DH - Global Total Least Squares misfit computation.
% MGTLS2DH(D,WL,WR,DH) = || sqrtm(WL) * (D - DH) * sqrtm(WR)||_F
%
% See GTLS for the possible formats of WL and WR.

[sd,N] = size(d);
if length(wl(:)) == sd % EWGTLS2
  wl = wl(:);  % column vector
  wl = sqrt(wl);
  wr = wr(:)'; % row vector
  wr = sqrt(wr);
  m  = norm(wl(:,ones(1,N)).*(d - dh).*wr(ones(sd,1),:),'fro');  
else  
  m = norm(chol(wl) * (d - dh) * chol(wr),'fro');
end