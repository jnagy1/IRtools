function Y = multiplyOneRegion( psfMatData, X );
%
%           Y = multiplyOneRegion( psfMatData, X );
%
%  This function computes the multiplication of a point spread
%  function matrix times an image vector; that is,
%               
%       Y = psfMatrix * X
%
%  Here we assume the dimension of the PSF is essentially the same
%  as the dimension of the image.
%
%  Input:
%   psfMatData  -  complex array containing the matrix data, usually
%                  computed from onePsfMatrix.m
%            X  -  array containing the image to which the psf matrix 
%                  is to be multiplied
%
%  Output:
%          Y  -  contains the result after multiplication.
%

%  J. Nagy  1/6/02

%
%  First we determine if any padding is needed.
%
Xpad = padarray(X, size(psfMatData) - size(X), 'post');

%
%  Now we perform the multiplications.
%  Note that this should work for 2D and 3D images.
%

Y = ifftn( psfMatData .* fftn( Xpad ) );

[nx, ny, nz] = size( X );
Y = real( Y(1:nx, 1:ny, 1:nz) );
