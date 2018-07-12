function [nregions, rsize] = partition_info(imsize, padsize)
%
%      [nregions, rsize] = partition_info(imsize, padsize);
%
%  Compute number of subretions and subregion size for overlap-add
%  and overlap-save convolution.
%
%  Input:
%         imsize  -  size of the image
% psfMatDataSize  -  size of A.matdata, where A is a psfMatrix object.
%
%  Output:
%       nregions  -  number of regions
%          rsize  -  region sizes
%

%  J. Nagy  1/7/02

%
%  First compute the basic sizes ...
%
psfSize = 2*padsize;
psfSize(psfSize == 0) = 1;
rsize = min( psfSize, imsize );
nregions = ceil(imsize ./ rsize);

%
%  This will avoid having one very small region at the end ...
%
% nregions = nregions - round( (rsize .* nregions - imsize) ./ rsize );
% rsize = ceil(imsize ./ nregions);

