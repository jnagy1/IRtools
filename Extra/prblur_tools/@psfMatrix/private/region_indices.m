function [Iidx, Jidx, Kidx] = region_indices( nregions, rsize )
%
%      [Iidx, Jidx, Kidx] = region_indices( nregions, rsize );
%
%  This function computes the starting and ending
%  row and column indices subregions of an
%  m-by-n array.  
%
%  Input:
%         nregions  -  vector whose length is the same a s the number
%                      of dimensions of the image.  This vector specifies
%                      the number of subregions in each dimension.
%         rsize     -  size of each subregion.
%
%  Output:
%         Iidx      -  array with two columns, giving region indices for the
%                      first image dimension.  The first column specifies
%                      the index of the beginning of each region.  The second
%                      column specifies the index of the end of each region.
%         Jidx      -  array with two columns, giving region indices for the 
%                      second image dimension.
%         Kidx      -  array with two columns, giving region indices for the
%                      third image dimension.  If there is not a third dimension,
%                      an empty array is returned.
%

%  J. Nagy  11/13/01


Iidx = [(0:nregions(1)-1)'*rsize(1)+1, (1:nregions(1))'*rsize(1)];

Jidx = [(0:nregions(2)-1)'*rsize(2)+1, (1:nregions(2))'*rsize(2)];

Kidx = [(0:nregions(3)-1)'*rsize(3)+1, (1:nregions(3))'*rsize(3)];


