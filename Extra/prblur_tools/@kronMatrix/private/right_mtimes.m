function y = right_mtimes(x, k)
%  Helper function for mtimes.m
%
%  
%  if x is a scalar, use the property that
%  kron(A, B) * x = kron(x * A, B)
%
%  if x is a vector such that x' = vec(X) 
%  use the property that
%  x * kron(A, B) = vec(B' * X * A)'
%

% Spring 2002, created by R. Wright

% 9/2002 updated by L. Perrone for incorporation into new
% kronMatrix class


disp('right')
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

elseif ( xsize(2) == A1size(1)*B1size(1) ) 
   y = zeros(xsize(1), A1size(2)*B1size(2));
   for j = 1:xsize(1) % iterate through the rows of x
    for i = 1:l %iterate through each kron product in the sum
      Ai = A{i};
      Bi = B{i};
      Aisize = size(Ai);
      Bisize = size(Bi);
      tmp = reshape(x(j,:), Bisize(1), Aisize(1));
      tmp = (Bi.')*tmp*Ai;
      y(j,:) = y(j,:) + ( reshape(tmp, 1, Aisize(2)*Bisize(2) ) );
    end % end the i loop
   end % end the j loop

elseif B1size(1) == xsize(1) & xsize(2) == A1size(1)
    % assume all Ai are same size, ditto for Bi
   y = zeros( B1size(2), A1size(2) );
   for i = 1:l
      Bi = B{i};
      tmp = (Bi.')*x*A{i};
      y = y + tmp;
   end
   %y=y';
else
  error('Dimension mismatch')
end
