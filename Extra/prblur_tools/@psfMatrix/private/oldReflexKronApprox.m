function [A, B] = reflexKronApprox( K, m, n )
%
%       [A, B] = reflexKronApprox( K, m, n );
%
% computes a kronecker sum approximation to the blurring matrix K
% that arises from the input PSF under reflexive boundary conditions.
% This approximation is done a-la Nagy-Ng-Perrone (see paper for details).
%
%  Input:
%         K - psfMatrix object
%         m - size of matrix A (assumed square)
%         n - size of matrix B (assumed square)
%
%  Output:
%      Matrices A and B such that K \aprrox A \otimes B.
%

%  L. Perrone, 4/28/02

%  Modifications:  
%  5/25/02, J. Nagy 
%           Cosmetic changes to incorporate into RestoreTools 
%
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
c = zeros(mp,1);
c(1) = mp;
c(2:2:end) = 1;
R = chol( toeplitz(c) );

Phat = R*PSF*R';

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
a = R \ ( U(:,1) * sqrt(S(1,1)) );
b = R \ ( V(:,1) * sqrt(S(1,1)) );

%
%  This construction corresponds to lexicographical ordering,
%  but kronMatrix does everything corresponding to vec ordering.
%  Thus, these two do not work together.
%%
%%A = build_toep(a, center(1), n) + buildHank(a, center(1), n);
%%B = build_toep(b, center(2), m) + buildHank(b, center(2), m);

A = build_toep(b, center(2), m) + buildHank(b, center(2), m);
B = build_toep(a, center(1), n) + buildHank(a, center(1), n);