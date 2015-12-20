function [m,dh] = mgtls2p(d,wl,wr,p)
% MGTLS2P - Global Total Least Squares misfit computation.
% [M,DH] = MGTLS2P(D,WL,WR,P) gives the GTLS2 misfit M and 
% the GTLS2 approximation DH of the data D by the model image(P).
%
% See GTLS for the possible formats of WL and WR.

% The computation is based on a transformation to a modified 
% TLS problem. As an alternative, see Corollary 2.13.

[sd,N] = size(d);

% Modify P and D
if length(wl(:)) == sd % EWGTLS2
  wl  = wl(:); % make it a column vector
  swl = sqrt(wl);
  wr  = wr(:)'; % make it a row vector
  swr = sqrt(wr);
  p   = swl(:,ones(1,size(r,1))) .* p;
  d   = swl(:,ones(1,N)) .* d .* swr(ones(sd,1),:); 
else
  swl = chol(wl);
  swr = chol(wr);
  p   = swl * p;
  d   = swl * d * swr;
end

% Compute the TLS misfit for the modified P and D
if nargout == 1
  m = mtlsp(d,p);
else
  [m,dh] = mtlsp(d,p);
end