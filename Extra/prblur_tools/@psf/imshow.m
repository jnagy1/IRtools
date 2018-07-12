function imshow( PSF )
%
%  Overload imshow for psf object.  This is preliminary
%  version, and needs some work.
%

% J. Nagy 1/8/02

I = PSF.image;
switch length(size(I{1}))
case 1
  plot(I)
case 2
  [m, n] = size(I);
  for i=1:m
    for j = 1:n
      subplot(m,n,(i-1)*n+j), imshow(I{i,j},[])
    end
  end
case 3
  [m, n, p] = size(I);
  if p > 1
    disp('We only display one of the spatially variant PSFs')
  end
  P = I{1,1,1};
  for k = 1:size(P,3)
    imshow(P(:,:,k),[min(P(:)), max(P(:))]), drawnow
  end
otherwise
  error('illegal PSF dimension')
end