function T = build_toep( c, k, n )
%
%  T = build_toep( c, k, n );
%
%  Given:
%    c - the nonzero part of a central column of a banded Toeplitz
%        matrix
%    k - index of c containing the diagonal element of T
%    n - dimension of the banded Toeplitz matrix
%
%  The banded Toeplitz matrix is constructed explicitly.
%

%  J. Nagy  2/11/02

m = length( c );

col = zeros(n,1);
row = col';
col(1:m-k+1,1) = c(k:m);
row(1,1:k) = c(k:-1:1)';
T = toeplitz( col, row );
