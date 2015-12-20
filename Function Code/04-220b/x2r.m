function r = x2r(x)
% X2R - Transformation from input/output to kernel representation.
% X2R(X) is a minimal kernel representation of B_io(X).

r = [x', -eye(size(x,2))];