function psfMatData = onePsfMatrix( PSF, center )
%
%       psfMatData = onePsfMatrix( PSF, center );
%
%  Construct matrix data for a single PSF.
%
%  Input:
%         PSF  -  array containing the PSF
%      center  -  array containing indices of center of the PSF.
%
%  Output:
%    psfMatData -  array containing the (complex) eigenvalues of the extended
%                  PSF
%

%  J. Nagy 1/6/02
%  3/22/07 got rid of PSF normalization.

%
%  Find the next power of 2 dimension of the extended PSF.
%

psf_size = size(PSF);
for i = 1:length(psf_size);
  e(i) = nextpow2( psf_size(i) );
end
m = 2 .^ (e+1);

%
%  Normalize the PSF so that the minimum pixel values is 0, and the 
%  sum of all entries is one (conservation of energy property).
%
%  Better to leave this to the user -- in the case of using
%  Tikhonov regularization, with regularization kernel 
%  [0 -1 0;-1 4 -1;0 -1 0]
%  you don't want to do this scaling.
%%%PSF = PSF - min(PSF(:));
%%%PSF = PSF / sum(PSF(:));

P = padarray(PSF, m-size(PSF), 'post');

P = circshift(P, -(center-1));

psfMatData = fftn(P);
