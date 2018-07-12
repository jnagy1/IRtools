function P = log( PSF )
%
%  Overload log for psf object.  This is preliminary
%  version, and needs some work.
%

%  J. Nagy 1/8/02

I = PSF.image;
c = PSF.center;
for k = 1:size(I,3)
  for i = 1:size(I,1)
    for j = 1:size(I,2)
      I{i,j,k} = log(I{i,j,k});
    end
  end
end
P = psf(I, c);