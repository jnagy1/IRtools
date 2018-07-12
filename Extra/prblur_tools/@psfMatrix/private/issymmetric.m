function YN = issymmetric(A)
%
%     issymmetric(A);
%
%  This function checks to see a PSF matrix is (close to) symmetric.  
%  It is used in svd.m to determine the best approach for computing SVD.
%
%  Returns 1 if A is symmetric, and 0 otherwise.
%
%

%  J. Nagy, 6/4/02

% Modifications:
% 11/12/02, J. Nagy
%           Now checks for circular symmetry.
%          
%


switch A.type
 
case 'variant'
  YN = 0;
case 'invariant'
  P = A.psf;
  c = P.center;
  c = c{1};
  P = P.image;
  P = P{1};
  n = size(P);

  P = padarray(P, max(n-2*c+1,0), 'pre');
  P = padarray(P, max(2*c-n-1,0), 'post');

  switch ndims(P)
  case 1
    Ps = (P + flipdim(P,1))/2;
  case 2
    Ps = (P + flipdim(P,1) + flipdim(P,2) + flipdim(flipdim(P,1),2))/4;
  case 3
    Ps = (P + flipdim(P,1) + flipdim(P,2) + flipdim(P,3) + ...
              flipdim(flipdim(P,1),2) + flipdim(flipdim(P,1),3) + ...
              flipdim(flipdim(P,2),3) + ...
              flipdim(flipdim(flipdim(P,1),2),3)) / 8;
  otherwise
    error('Can only do 1, 2 or 3 dimensions.')
  end

  if ( norm(Ps(:) - P(:)) <= sqrt(eps)*norm(P(:)) )
    YN = 1;
  else
    YN = 0;
  end

otherwise
  error('Invalid psfMatrix type')
end
