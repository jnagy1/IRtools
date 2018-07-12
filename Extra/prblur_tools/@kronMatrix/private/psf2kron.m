function varargout = kronApprox2( varargin )
%
%  Compute a Kronecker product approximation of a psfMatrix, K:
%           K \approx A \otimes B
%
%  On Entry: 
%           PSF
%           g - the right hand side vector (blurred image)
%
%  On Exit:
%      * If there is only one ouput argument, then a kronMatrix object
%        is returned. (See help kronMatrix for more inforation.)
%      * If there are two output arguments, then the matrices A and B
%        are returned.
%

PSF = varargin{1};
if nargin == 1
  imgDimension = size(PSF);
  mm = imgDimension(1);
  nn = imgDimension(2);
else
  imgDimension = size(varargin{2});
  mm = imgDimension(1);, nn = imgDimension(2);
end
if length( imgDimension ) ~= 2 
  error('kronApprox works only for 2D problems')
end

[c1, c2] = find(PSF == max(PSF(:)));
center = [c1, c2];

[m,n] = size(PSF);
if ( m ~= n )
  error('For now, we expect PSF to be square')
end

%
% Compute weighted PSF.
%
for i = 1:center(1)
  Aweights(i,1) = sqrt( i+m-center(1) );
end;
for i = center(1)+1:m
  Aweights(i,1) = sqrt( m+center(1)-i);
end;
for i = 1:center(2)
  Bweights(i,1) = sqrt( i+n-center(2) );
end;
for i = center(2)+1:n
  Bweights(i,1) = sqrt( n+center(2)-i );
end;
Phat = (Aweights*Bweights').*PSF;
%
% Compute SVD of weighted PSF, which is then used to construct
% the separable approximation.
%
[U,S,V] = svd( Phat );


%
% check to make sure first column looks like
% a Gaussian, and is not inverted.
%
minU = abs(min(min(U(:,1))));
maxU = max(max(abs(U(:,1))));
if minU == maxU
  U = -U;
  V = -V;
end

%
% Construct approximation.
%
M(:,1) = ( U(:,1) * sqrt(S(1,1)) ) ./ Aweights;
M(:,2) = ( V(:,1) * sqrt(S(1,1)) ) ./ Bweights;

A = build_toep( M(:,1), c1, mm );
B = build_toep( M(:,2), c2, nn );

if nargout == 1
  varargout{1} = kronMatrix(A, B);
elseif nargout == 2
  varargout{1} = A;
  varargout{2} = B;
else
  error('Incorrect number of output arguments')
end  