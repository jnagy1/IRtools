function P = adjustPsfSize(P)
%
%      P = adjustPsfSize(P);
%
%  If the size of a PSF is small compared to the dimensions of
%  the image, then our matrix vector multiply functions will be
%  inefficient due to our overlap add and overlap save approach.
%  Therefore, we make the minimum PSF size to be 32.
%
%  Input:  P - psf object
%  Output: P - psf object
%

%  J. Nagy 1/12/02

I = P.image;
center = P.center;

for k = 1:size(I,3)
  for j = 1:size(I,2)
    for i = 1:size(I,1)
      if length(I{i,j,k}) < 32
        I1 = zeros(32*ones(1,length(size(I{i,j,k}))));
        I1(1:size(I{i,j,k},1), 1:size(I{i,j,k},2), 1:size(I{i,j,k},3)) = I{i,j,k};
        I{i,j,k} = I1;
      end
    end
  end
end
P = psf(I, center);