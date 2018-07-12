%
%  T = build_toep( c, d, n );
%
%  Given:
%    c - the nonzero part of a central column of a banded Toeplitz
%        matrix
%    d - index of c containing the diagonal element of T
%    n - dimension of the banded Toeplitz matrix
%
%  The banded Toeplitz matrix is constructed explicitly.
%
function T = build_toep( c, d, n )

m = length( c );

col = zeros(n,1);
row = col';
col(1:m-d+1,1) = c(d:m);
row(1,1:d) = c(d:-1:1)';
T = toeplitz( col, row );
