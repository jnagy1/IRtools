function [Iidx, Jidx, Kidx] = eregion_indices(RIidx, RJidx, RKidx, rsize )
%
%      [Iidx, Jidx, Kidx] = eregion_indices( RIidx, RJidx, RKidx, rsize );
%
%  This function computes the starting and ending
%  row and column indices subregions of an
%  m-by-n array.  
%
%  Input:
%         RIidx   \  
%         RJidx     -  indices coming from region_indices.m
%         RKidx   /
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

%  J. Nagy 11/27/01

%
%  Assume we have at least one dimension
%
[m, n] = size(RIidx);
Iidx = RIidx + [zeros(m,1), 2*fix(rsize(1)/2)*ones(m,1)];

%
% If there is a second dimension, ...
%
if rsize(2) > 1
  [m, n] = size(RJidx);
  Jidx = RJidx + [zeros(m,1), 2*fix(rsize(2)/2)*ones(m,1)];
else
  Jidx = [1,1];
end

%
% Finally, if there is a third dimension, ...
%
if rsize(3) > 1
  [m, n] = size(RKidx);
  Kidx = RKidx + [zeros(m,1), 2*fix(rsize(3)/2)*ones(m,1)];
else
  Kidx = [1,1];
end

