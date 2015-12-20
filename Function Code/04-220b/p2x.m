function x = p2x(p,tol)
% P2X - From image to an input/output representation.
% P2X(P) is an input/output representation of image(P).
% P2X(P,TOL) ill conditioning of the transformation is 
% checked with the user defined tolerance TOL.

if nargin == 1
  tol = 1e-14; % default 
end
  
% Find a minimal image representation
p = minp(p); 

% Define the number of inputs
[sd,m] = size(p);

% Check the conditioning of the transformation
if rcond(p(1:m,:)) < tol
  error('Conversion to an I/O representation is ill conditioned.')
end

% Transform
x = ( p(m+1:end,:) / p(1:m,:) )';