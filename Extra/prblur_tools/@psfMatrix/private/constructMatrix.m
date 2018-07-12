function psfMatData = constructMatrix( PSFs, center )
%
%       psfMatData = constructMatrix( PSFs, center );
%
%  Construct psfMatrix data.
%
%  Given several PSFs and the locations of the corresponding point sources,
%  this function sets up the data needed to do efficient matrix-vector
%  multiplication (convolution) with a possibly space variant PSF.
%
%  Input:
%        PSFs  -  d-dimensional cell array containing the PSFs
%                 e.g., for 2-d problems, it is a 2-d cell array,
%                       for 3-d problems, it is a 3-d cell array.
%      center  -  cell array having same shape as PSFs, containing
%                 location of the center of the corresponding PSF.
%
%  Output:
%   psfMatData -  cell array containing the (complex) data needed to
%                 do efficient matrix-vector multiplications.
%

%  J. Nagy  1/12/02

psfMatData = cell(size(PSFs));

for k = 1:size(PSFs,3)
  for i = 1:size(PSFs,1)
    for j = 1:size(PSFs,2)

       psfMatData{i,j,k} = onePsfMatrix(PSFs{i,j,k}, center{i,j,k});

    end
  end
end
