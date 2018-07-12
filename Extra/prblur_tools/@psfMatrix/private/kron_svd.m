function [U, s, V] = kron_svd(A)
%
%   [U, s, V] = kron_svd();
%
%  Given a kronMatrix object, this function computes
%  an "SVD" factorization that is to be used for a
%  psfMatrix object.
%
%  On Entry:
%     A  -  kronMatrix obect
%
%  On Exit:
%     U, V  -  kronMatrix objects, orthogonal
%        s  -  column vector containing the singular values
%              of A.  Note that these are not sorted from largest
%              smallest.
%
[U, S, V] = svd(A);

if isa(S.a, 'cell')
  s = kron(diag(S.a{1}), diag(S.b{1}));
  for i = 2:length(S.a)
    s = s + kron(diag(S.a{i}), diag(S.b{i}));
  end
else
  s = kron(diag(S.a), diag(S.b));
end
