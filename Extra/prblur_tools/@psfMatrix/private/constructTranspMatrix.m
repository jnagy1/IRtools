function [psfTranspMatData] = constructTranspMatrix( PSFs, center )
%
%       psfTranspMatData = constructTranspMatrix( PSFs, center );
%
%  Construct tranpose psfMatrix data.
%
%  Given several PSFs and the locations of the corresponding point sources,
%  this function sets up the data needed to do efficient matrix-vector
%  multiplication with the transpose with a possibly space variant PSF.
%
%  Input:
%        PSFs  -  d-dimensional cell array containing the PSFs
%                 e.g., for 2-d problems, it is a 2-d cell array,
%                       for 3-d problems, it is a 3-d cell array.
%      center  -  cell array having same shape as PSFs, containing
%                 location of the center of the corresponding PSF.
%
%  Output:
%   psfTranspMatData -  cell array containing the (complex) data needed to
%                 do efficient matrix-vector multiplications with transpose.
%

%  J. Nagy  7/2/2012

psfTranspMatData = cell(size(PSFs));

for k = 1:size(PSFs,3)
  for i = 1:size(PSFs,1)
    for j = 1:size(PSFs,2)
       PSFrot = flipdim(flipdim(flipdim(PSFs{i,j,k},1),2),3);
       center_rot = size(PSFrot) - center{i,j,k} + 1;
       psfTranspMatData{i,j,k} = onePsfMatrix(PSFrot, center_rot);

    end
  end
end
