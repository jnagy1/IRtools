function [TIidx, TJidx, TKidx] = tregion_indices(RIidx, RJidx, RKidx, padsize)
%
%    [TIidx, TJidx, TKidx] = tregion_indices(RIidx, RJidx, RKidx, padsize);
%
%  This is needed in variantTransposeMultiply.m
%
%  Input:
%         RIidx   \
%         RJidx     -  indices coming from region_indices.m
%         RKidx   /
%       padsize     -  amount of padding for boundary conditions
%
%  Output:
%        TIidx      -  array with two columns, giving region indices for the
%                      first image dimension.  The first column specifies
%                      the index of the beginning of each region.  The second
%                      column specifies the index of the end of each region.
%        TJidx      -  array with two columns, giving region indices for the
%                      second image dimension.
%        TKidx      -  array with two columns, giving region indices for the
%                      third image dimension. 
%

%  J. Nagy 1/9/02


TIidx = ones(size(RIidx));
TIidx(2:end,1) = padsize(1)+1;
TIidx(1:end-1,2) = RIidx(1,2);
TIidx(end,2) = RIidx(1,2) + padsize(1);

TJidx = ones(size(RJidx));
TJidx(2:end,1) = padsize(2)+1;
TJidx(1:end-1,2) = RJidx(1,2);
TJidx(end,2) = RJidx(1,2) + padsize(2);

TKidx = ones(size(RKidx));
TKidx(2:end,1) = padsize(3)+1;
TKidx(1:end-1,2) = RKidx(1,2);
TKidx(end,2) = RKidx(1,2) + padsize(3);