function [U, s, V] = dct_svd(PSF, center)
%
%   [U, s, V] = dct_svd(PSF, center);
%
%  Given a symmetric PSF this function computes an "SVD" factorization
%  for reflexive boundary conditions.  That is, the "SVD" is the
%  "Spectral Value Decomposition":
%            A = USV',
%  where U = V = inverse DCT matrix.
%
%  On Entry:
%     A  -  psfMatrix obect, with periodic boundary conditions
%
%  On Exit:
%     U, V  -  transformMatrix objects, with 
%              U.transform = V.transform = 'dct'
%        s  -  column vector containing the eigenvalues
%              of A.  Note that these are not sorted from largest
%              smallest.
%
e1 = zeros(size(PSF));, e1(1,1) = 1;
c = center;, n = size(PSF);
if ( ndims(PSF) == 1 )
  s = dct(PSF) ./ dct(e1);
elseif ( ndims(PSF) == 2 )
  P = zeros(size(PSF));
  P(1:n(1)-c(1)+1,1:n(2)-c(2)+1) = PSF(c(1):n,     c(2):n       );
  P(1:n(1)-c(1)+1,1:n(2)-c(2)  ) = P(1:n(1)-c(1)+1,1:n(2)-c(2)  ) + PSF(c(1):n,  c(2)+1:n);
  P(1:n(1)-c(1),  1:n(2)-c(2)+1) = P(1:n(1)-c(1),  1:n(2)-c(2)+1) + PSF(c(1)+1:n,c(2):n  );
  P(1:n(1)-c(1),  1:n(2)-c(2)  ) = P(1:n(1)-c(1),  1:n(2)-c(2)  ) + PSF(c(1)+1:n,c(2)+1:n);
  s = dct2(P) ./ dct2(e1);
  s = s(:);
else
  error('This does not work for more than two dimensions.')
end

U = transformMatrix('dct');
U = U';
V = transformMatrix('dct');
V = V';
