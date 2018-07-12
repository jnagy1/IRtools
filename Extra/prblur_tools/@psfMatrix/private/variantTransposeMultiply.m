function Y = variantTransposeMultiply( psfMatData, X, padsize );
%
%           Y = variantTransposeMultiply( psfMatData, X, padsize );
%
%  This function computes the transpose multiplication of a spatially variant
%  point spread function (PSF) times an image:
%                y = A'*x
%
%  Here we assume A is made up of a several space invariant PSFs, and
%  use piece-wise constant interpolation of them to define the spatially
%  variant PSF.  That is, A has the form:
%
%       A  = D1*A1  + D2*A2  + ... + Dp*Ap
%  ==>  A' = A1'*D1 + A2'*D2 + ... + Ap'*Dp
%
%  Input:
%   psfMatData  -  cell array containing the matrix data of each of the
%                  individual PSFs.  This matrix data is usuall computed
%                  from onePsfMatrix.m
%            X  -  array containing the image to which the psfMatrix
%                  is to be multiplied.
%
%  Output:
%            Y  -  contains the result after PSF transpose multiplication.
%

%  J. Nagy  1/7/02

imsize = size( X ) - 2*padsize;

%
%  We partition the image domain into regions of equal sizes,
%  according to the number of PSFs we have ...
%
nregions = size(psfMatData);
rsize = ceil(imsize ./ nregions);

%
%  In order for this to be consistent for 2-D and 3-D images, we need to make
%  sure there is a third dimension ...
%
if length(imsize) == 1
  imsize = [imsize, 1, 1];
  rsize = [rsize, 1, 1];
  nregions = [nregions, 1, 1];
  padsize = [padsize, 0, 0];
elseif length(imsize) == 2
  imsize = [imsize, 1];
  rsize = [rsize, 1];
  nregions = [nregions, 1];
  padsize = [padsize, 0];
end

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

%
%  Now we get information about beginning and ending indices of subregions
%  so we can "put" and "get" subregions correctly ...
%
[RIidx, RJidx, RKidx] = region_indices( nregions, rsize );
[EIidx, EJidx, EKidx] = eregion_indices( RIidx, RJidx, RKidx, 2*padsize );
RIidx = RIidx + padsize(1);, RIidx(1) = 1;, RIidx(end) = RIidx(end) + padsize(1);
RJidx = RJidx + padsize(2);, RJidx(1) = 1;, RJidx(end) = RJidx(end) + padsize(2);
RKidx = RKidx + padsize(3);, RKidx(1) = 1;, RKidx(end) = RKidx(end) + padsize(3);
[TIidx, TJidx, TKidx] = tregion_indices(RIidx, RJidx, RKidx, padsize);

%
%  Now loop over all the subregions ...
%

Y = zeros( size(X) );
for k = 1:nregions(3)
  for j = 1:nregions(2)
    for i = 1:nregions(1)

      eregionSize = [diff(EIidx(i,:)), diff(EJidx(j,:)), diff(EKidx(k,:))] + 1;
      Xt = zeros( eregionSize );

      Xt(TIidx(i,1):TIidx(i,2), TJidx(j,1):TJidx(j,2), TKidx(k,1):TKidx(k,2)) =  ...
            X(RIidx(i,1):RIidx(i,2), RJidx(j,1):RJidx(j,2), RKidx(k,1):RKidx(k,2));

      p = padsize(1:length(size(Xt)));
      Xt = padarray(Xt, p,'both');

      Y(EIidx(i,1):EIidx(i,2), EJidx(j,1):EJidx(j,2), EKidx(k,1):EKidx(k,2)) = ...
            Y(EIidx(i,1):EIidx(i,2), EJidx(j,1):EJidx(j,2), EKidx(k,1):EKidx(k,2)) + ...
            invariantMultiply( psfMatData{i,j,k}, Xt, p );
            %invariantMultiply( conj(psfMatData{i,j,k}), Xt, p );
 
    end
  end
end
idx = [padsize+1; padsize+imsize];
Y = Y( idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3) );

