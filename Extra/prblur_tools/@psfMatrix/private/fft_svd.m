function [U, s, V] = fft_svd(PSF, center)
%
%   [U, s, V] = fft_svd(PSF, center);
%
%  Given a PSF this function computes an "SVD" factorization
%  for periodic boundary conditions.  That is, the "SVD" is
%  "Spectral Value Decomposition":
%            A = USV',
%  where U = V = inverse DFT matrix.
%
%  On Entry:
%     A  -  psfMatrix obect, with periodic boundary conditions
%
%  On Exit:
%     U, V  -  transformMatrix objects, with 
%              U.transform = V.transform = 'fft'
%        s  -  column vector containing the eigenvalues
%              of A.  Note that these are not sorted from largest
%              smallest.
%
P = circshift(PSF, -(center - 1));
S = fftn(P);
s = S(:);

U = transformMatrix('fft');
U = U';
V = transformMatrix('fft');
V = V';