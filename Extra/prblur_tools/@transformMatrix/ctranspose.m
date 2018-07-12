function A = ctranspose(A);
%
%  CTRANSPOSE The transpose of the trasformMatrix matrix 
%

%  J. Nagy 6/2/02

if A.transpose == 0
  A.transpose = 1;
else
  A.transpose = 0;
end
