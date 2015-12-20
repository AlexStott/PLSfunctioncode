function m = mtlsdh(d,dh)
% MTLSDH - Total Least Squares misfit computation.
% MTLSDH(D,DH) = || D - DH ||_F

m = norm(d-dh,'fro');