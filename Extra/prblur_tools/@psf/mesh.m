function mesh( PSF )
%
%  Overload mesh for psf object.  This is preliminary
%  version, and needs some work.
%

%  J. Nagy 1/8/02

I = PSF.image;
switch length(size(I{1}))
case 1
  plot(I)
case 2
  [m, n] = size(I);
  for i=1:m
    for j = 1:n
      subplot(m,n,(i-1)*n+j), mesh(I{i,j})
    end
  end
case 3
  if size(I,3) > 1
    disp('We only display one of the spatially variant PSFs')
  end
  P = I{1,1,1};
  v = [0, size(P,1), 0, size(P,2), min(P(:)), max(P(:))];
  for k = 1:size(P,3)
    mesh(P(:,:,k)), axis(v), drawnow
  end
otherwise
  error('illegal PSF dimension')
end