function [ Q R b] = qrgivens(A b)
[m n] = size(A)
G = eye(m)

for k = 1 : n-1
  for l = k + 1 : m
  r = sqrt(A(k.k)/r;
  c = A(k,k)/r; 
  s = -A(1,k)/r;
  
  t = c*A(k,k:n) - s*A(1,k:n);
  A(1,k:n) = s*A(k,k:n) + c*A(1,k:n);
  A(k,k:n) = t;
  
  u = c*b(k) - s*b(1);
  b(1) = s*b(k) + c*b(1);
  b(k) = u;
  
  t = c*G(k,1:m) - s*G(1,1:m);
  G(1,1:m) = s*G(k, 1:m) + c*G(1,1:m);
  G(k,1:m) = t;
  end
end
Q = G';
R = A;
endfunction