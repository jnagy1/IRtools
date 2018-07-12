function Y = invariantMultiply( psfMatData, X, padsize );
%
%           Y = invariantMultiply( psfMatData, X, padsize );
%
%  This function computes the multiplication of a spatially invariant
%  point spread function (PSF) times an image (i.e., convolution):
%                y = A*x
%
%  Here we assume A is made up of a single PSF whose extent may be
%  (much) smaller than the image.
%
%  Input:
%   psfMatData  -  complex array containing the matrix data, usually
%                  computed from onePsfMatrix.m
%            X  -  array containing the image to which the psfMatrix
%                  is to be multiplied.
%
%  Output:
%            Y  -  contains the result after PSF multiplication.
%

%  J. Nagy  1/7/02

imsize = size( X ) - 2*padsize;

%
%  In order for this to be consistent for 2-D and 3-D images, we need to make 
%  sure there is a third dimension ...
%
if length(imsize) == 1
  imsize = [imsize, 1, 1];
  padsize = [padsize, 0, 0];
elseif length(imsize) == 2
  imsize = [imsize, 1];
  padsize = [padsize, 0];
end
 
%
%  partition_info computes number of subregions, and their sizes
%
[nregions, rsize] = partition_info(imsize, padsize);

%
%  Coding the rest of this will be easier if all of the image subregions
%  have the same dimensions.  If it's not, we pad with a few zeros to make
%  it so ...
%
padsize1 = rsize .* nregions - imsize;
if any( padsize1 < 0 )
  error('Something is wrong here ...')
end
X = padarray(X, padsize1, 'post');

Y = zeros( size(X) );

%
%  Now we get information about beginning and ending indices of subregions
%  so we can "put" and "get" subregions correctly ...
%
[RIidx, RJidx, RKidx] = region_indices( nregions, rsize );
[EIidx, EJidx, EKidx] = eregion_indices( RIidx, RJidx, RKidx, 2*padsize );
Tidx = [padsize+1; padsize+rsize];

%
%  Now loop over all the subregions ...
%
for k = 1:nregions(3)
  for j = 1:nregions(2)
    for i = 1:nregions(1)
      Xt = X(EIidx(i,1):EIidx(i,2), EJidx(j,1):EJidx(j,2), EKidx(k,1):EKidx(k,2));

      Yt = multiplyOneRegion( psfMatData, Xt );
      
      Y(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2)) = ...
           Yt(Tidx(1,1):Tidx(2,1), Tidx(1,2):Tidx(2,2), Tidx(1,3):Tidx(2,3));
    end
  end
end


Y = Y(1:imsize(1), 1:imsize(2), 1:imsize(3));


