function A = ctranspose(A);
%
%  CTRANSPOSE The transpose of the psfMatrix matrix is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  J. Nagy 5/2/01

if A.transpose == 0
  A.transpose = 1;
else 
  A.transpose = 0;
end