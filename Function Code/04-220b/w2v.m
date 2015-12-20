% W2V - Convert weight matrices to covariance matrices.
% V = W2V(W,C)
% C == 1 (EWTLS): V := 1 ./ W
% C == 2 (WTLS) : V(:,:,i) := inv(W(:,:,i)), i = 1,...,N
% C == 3 (FTLS) : V := inv(W)

function w = w2v(w,c)

if c == 1
  w = 1./w;
elseif c == 2
  N = size(w,3);
  for i = 1:N
    w(:,:,i) = inv(w(:,:,i));
  end
else
  w = inv(w);
end
