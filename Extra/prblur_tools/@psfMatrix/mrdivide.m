function A = mrdivide(A,arg2)
%
%   Overload Matrix scalar division  for psfMatrix
%
%   Implements A/s where A is psfMatrix and s is a scalar
%
%   Result returned is a psfMatrix
%

%  J. Nagy & K. Lee  1/18/02

M1=A.matdata;
for k=1:size(M1,3)
  for j=1:size(M1,2)
    for i=1:size(M1,1)
      M1{i,j,k}=M1{i,j,k}/arg2;
    end
  end
end
A.matdata=M1;