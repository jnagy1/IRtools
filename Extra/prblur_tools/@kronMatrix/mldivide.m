function x = mldivide(K,y)
%  kronMatrix/mldivide
%
%     solves the linear system of equations
%     Kx = y using only the first kron product in the
%     sum for K.  Returns a double.
%

% R. Wright Spring 2002
% L. Perrone 10/2002 updated for new kronMatrix class


A = K.a{1};
B = K.b{1};
[ma, na] = size(A);
[mb, nb] = size(B);
Y = reshape(y, mb, ma);
Z = A \ Y';
X = B \ Z';
x = reshape(X, nb*na, 1);

if length(K) > 1
  disp('A loss of accuracy may have resulted; only one kronecker product was used to perform the back-solve.')
end