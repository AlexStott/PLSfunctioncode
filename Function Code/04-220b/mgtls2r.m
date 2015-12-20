function [m,dh] = mgtls2r(d,wl,wr,r)
% MGTLS2R - Global Total Least Squares misfit computation.
% [M,DH] = MGTLS2R(D,WL,WR,R) gives the GTLS2 misfit M and 
% the GTLS2 approximation DH of the data D by the model ker(R).
%
% See GTLS for the possible formats of WL and WR.

% The computation is based on a transformation to a modified 
% TLS problem. As an alternative, see Corollary 2.11.

[sd,N] = size(d);

% Modify R and D
if length(wl(:)) == sd % EWGTLS2
  wl  = wl(:); % make it a column vector
  swl = sqrt(wl);
  wr  = wr(:)'; % make it a row vector
  swr = sqrt(wr);
  r   = r ./ swl(:,ones(1,size(r,1)))';
  d   = swl(:,ones(1,N)) .* d .* swr(ones(sd,1),:); 
else 
  swl = chol(wl);
  swr = chol(wr);
  r   = r / swl;
  d   = swl * d * swr;
end

% Compute the TLS misfit for the modified R and D
if nargout == 1
  m = mtlsr(d,r);
else
  [m,dh] = mtlsr(d,r);
end