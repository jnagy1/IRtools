function y = left_mtimes(k, x)

%  y = left_mtimes(k, x)
%  Helper function for mtimes.m
%  k is a kronMatrix object, x is a double (scalar,vector,or matrix)
%
%  if x is a scalar, use the property that
%  x * kron(A, B) = kron(x * A, B)
%
%  if x is a vector such that x = vec(X),
%  use the property that 
%  kron(A, B) * x = vec(B * X * A')
%

% Spring 2002, created by R. Wright

% 9/2002 updated by L. Perrone for incorporation into new
% kronMatrix class


A = k.a;
B = k.b;
l = length(A);
xsize = size(x);
A1size = size(A{1});
B1size = size(B{1});

if length(x) == 1
   for i=1:l
      A{i} = A{i}*x;
   end
   y = kronMatrix(A,B);

elseif xsize(1) == A1size(2)*B1size(2)
   y = zeros(A1size(1)*B1size(1) , xsize(2) );
   for j = 1:xsize(2) % iterate through the columns of x
    for i = 1:l % iterate through each kron product in the sum
      Ai = A{i};
      Bi = B{i};
      Aisize = size(Ai);
      Bisize = size(Bi);
      tmp = reshape(x(:,j), Bisize(2), Aisize(2));
      tmp = Bi*tmp*(Ai.');
      y(:,j) = y(:,j) + reshape(tmp, Aisize(2)*Bisize(2), 1);
    end % end the i loop
   end % end the j loop

elseif B1size(2) == xsize(1) & xsize(2) == A1size(2)
    % assume all Ai are same size, ditto for Bi
   y = zeros( B1size(1), A1size(1) );
   for i = 1:l
      Ai = A{i};
      tmp = B{i}* x * Ai.';
      y = y + tmp;  
   end
else
  error('Dimension mismatch')
end
