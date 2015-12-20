function p = x2p(x)
% X2P - Transformation from input/output to image representation.
% X2P(X) is a minimal image representation of B_io(X).

p = [eye(size(x,1)); x'];