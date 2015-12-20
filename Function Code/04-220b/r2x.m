function x = r2x(r,tol)
% R2X - From kernel to an input/output representation.
% R2X(R) is an input/output representation of ker(R).
% R2X(R,TOL) ill conditioning of the transformation is 
% checked with the user defined tolerance TOL.

if nargin == 1
  tol = 1e-14; % default 
end

% Find a minimal kernel representation
r = minr(r); 

% Define number of inputs and outputs
[p,sd] = size(r);
m = sd - p;

% Check the conditioning of the transformation
if rcond(r(:,m+1:end)) < tol
  error('Conversion to an I/O representation is ill conditioned.')
end

% Transform
x = -( r(:,m+1:end) \ r(:,1:m) )';