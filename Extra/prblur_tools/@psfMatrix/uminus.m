function A  = uminus( A )
%
%  Overload uminus for a psfMatrix object.
%  Given A, this retruns -A.
%

%  K. Lee  1/9/02

M = A.matdata;
for k = 1:size(M,3);
  for j = 1:size(M,2);
    for i = 1:size(M,1);
      M{i,j,k} = -M{i,j,k};
    end
  end
end
A.matdata = M;
