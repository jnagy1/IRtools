function [m,n]=imageSize(A)
%
%  Tells the total image size of the psfMatrix object
%
%

%  K. Lee 1/28/02

A=A.psf;
A=A.image;
[rows,cols]=size(A);
m=0;
n=0;
for i = 1:cols
  A1=A{i};
  [r,A1cols]=size(A1);
  n=n+A1cols;
end
for i = 1:rows
  A1=A{i};
  [r,A1cols]=size(A1);
  m = m+r;
end