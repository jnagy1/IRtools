function display(P)
%
%  DISPLAY(P) Display information about the psf object.
%

%  J. Nagy 5/1/01

stg_img = sprintf('There are %d PSFs', prod(size(P.image)));
disp(stg_img)