function G = kronApprox(K, g, Rtol)
%
%  Compute a Kronecker product approximation of a psfMatrix, K:
%           K \approx G = \sum [ A{i} (x) B{i} ]
%
%  On Entry: 
%           K - psfMatrix object
%           g - the right hand side vector (blurred image)
%        Rtol - truncation tolerance for the number of terms in the
%               sum of Kronecker product approximation
%                OR
%                positive integer specifying the number of terms in the
%                approximation.
%
%  On Exit:
%      G -  kronMatrix object (See help kronMatrix for more inforation.)
%

%  L. Perrone, 5/2003

% Modifications:
% 6/2003, J. Nagy
%         Combined some codes from original kronMatrix.m with
%         original kronApprox.
%

switch nargin
case 1
  g = [];
  Rtol = [];
case 2
  Rtol = [];
end

switch K.type

case 'invariant'

  P1 = K.psf;
  P2 = P1.image;
  PSF = P2{1};
  psfDimension = size(PSF);

  if isempty(g)
    m = psfDimension(1);
    n = psfDimension(2);
  else
    imgDimension = size(g);
    m = imgDimension(1);
    n = imgDimension(2);
  end
  if length( psfDimension ) ~= 2 
    error('kronApprox works only for 2D problems')
  end

  if isempty(Rtol)
    Rtol = 0.01;
  end
    
  switch K.boundary
  case 'zero'
    [A, B] = zeroKronApprox(K, m, n, Rtol);
  case 'reflexive'
    [A, B] = reflexKronApprox(K, m, n, Rtol);
  case 'neumann'
    [A, B] = reflexKronApprox(K, m, n, Rtol);
  otherwise
    error('kronApprox not defined for given boundary condition')
  end
  G = kronMatrix(A, B);
  %if iscell(G.a)
  %  disp(sprintf('********** number of terms in Kronecker sum approximation is %d', length(G.a)))
  %else
  %  disp('************* number of terms in Kronecker sum approximation is 1')
  %end

case 'variant'

  if isempty(g)
    error('you must input blurred image for space variant blurs')
  else
    imgDimension = size(g);
    m = imgDimension(1);
    n = imgDimension(2);
  end
  if isempty(Rtol)
    Rtol = 0.01;
  end
  [A, B] = variantKronApprox(K, m, n, Rtol);
  G = kronMatrix(A, B);

otherwise

  error('Invalid psfMatrix type')

end

