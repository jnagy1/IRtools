function s = size( P )
%
%  Overload size for psf object.
%  Returns the size of the PSF image.
%  For spatially variant blurs, returns the max size
%  of all the PSFs.
%

%  J. Nagy  1/6/02

I = P.image;

if prod(size(I)) > 1;
  s = [];
  for k = 1:size(I,3)
    for i = 1:size(I,1)
      for j = 1:size(I,2)
        s = [s; size(I{i,j,k})];
      end
    end
  end
  s = max(s);
else
  s = size(I{1});
end
