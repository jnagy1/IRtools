function [A, B] = zeroKronApprox( K, m, n )
%
%       [A, B] = zeroKronApprox( K, m, n );
%
% computes a kronecker sum approximation to the blurring matrix K
% that arises from the input PSF under zero boundary conditions.
% This approximation is done a-la Kamm-Nagy (see paper for details).  
%
%  Input:
%         K - psfMatrix object
%         m - size of matrix A (assumed square)
%         n - size of matrix B (assumed square)
%
%  Output:
%      Matrices A and B such that K \aprrox A \otimes B.
%

%  J. Nagy  2/11/02

%  Modifications:  This used to be zeroSepMatrix, which did not compute
%                  A and B. 
%  5/25/02, J. Nagy 
%           This now computes A and B.  
%
%  11/17/02, J. Nagy
%            This was designed for image processing problems, where
%            it is common to use lexicographical (row) ordering.
%            But @kronMatrix functions where designed using vec (column)
%            ordering.  This inconsistency has been fixed.
 
P1 = K.psf;
P2 = P1.image;
PSF = P2{1};
c1 = P1.center;
center = c1{1};

[mp, np] = size(PSF);

if ( mp ~= np )
  error('For now, we expect PSF to be square')
end

%
% Compute weighted PSF.
%
for i = 1:center(1)
  Aweights(i,1) = sqrt( i+mp-center(1) );
end;
for i = center(1)+1:mp
  Aweights(i,1) = sqrt( mp+center(1)-i);
end;
for i = 1:center(2)
  Bweights(i,1) = sqrt( i+np-center(2) );
end;
for i = center(2)+1:np
  Bweights(i,1) = sqrt( np+center(2)-i );
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
a = ( U(:,1) * sqrt(S(1,1)) ) ./ Aweights;
b = ( V(:,1) * sqrt(S(1,1)) ) ./ Bweights;

%
%  This construction corresponds to lexicographical ordering,
%  but kronMatrix does everything corresponding to vec ordering.
%  Thus, these two do not work together.
%%
%%A = build_toep(a, center(1), n);
%%B = build_toep(b, center(2), m);

%  To fix this, we just need to switch A and B.  Then kronMatrix
%  can continue to be used corresponding to a vec ordering.
%
A = build_toep(b, center(2), m);
B = build_toep(a, center(1), n);

