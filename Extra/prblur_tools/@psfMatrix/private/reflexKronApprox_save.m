function [A, B] = reflexKronApprox( K, m, n, Rtol )
%
%       [A, B] = reflexKronApprox( K, m, n, Rtol );
%
% computes a kronecker sum approximation to the blurring matrix K
% that arises from the input PSF under reflexive boundary conditions.
% This approximation is done a-la Nagy-Ng-Perrone (see paper for details).
%
%  Input:
%         K - psfMatrix object
%         m - size of matrix A (assumed square)
%         n - size of matrix B (assumed square)
%         Rtol - relative tolerance ( <=1 ); determines number of kron products
%                in the approximation.  If si/s1 > Rtol, 
%                where si is the ith largest
%                singular value of the weighted PSF,
%                then Ai(x)Bi will be included in the approximation.
%                IF ONLY ONE TERM IN THE SUM IS DESIRED, USE Rtol = 1.
%
%  Output:
%      Cells A and B such that K \approx \sum_i[ A{i} \otimes B{i} ].
%

%  L. Perrone, 4/28/02

%  Modifications:  
%  5/25/02, J. Nagy 
%           Cosmetic changes to incorporate into RestoreTools 
%
%  9/??/02  L. Perrone
%           This code now conforms to the new kronMatrix class,
%           where fields K.a and K.b are cell arrays containing
%           (possibly) more than one matrix, depending on Rtol.
%  11/22/02 L. Perrone
%           The kronMatrix class should work for image processing
%           problems, where it is common to use lexicographical (row)
%           ordering.  However, the kronMatrix class was designed
%           using vec(column) ordering (unlike all other classes
%           in the RestoreTools package).  Until now, the inconsistency
%           remained hidden because the PSFs in use
%           had been symmetric or close to symmetric.  
%           This discrepancy has now been fixed.



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
A = cell(1);
B = cell(1);
i=1;
PP = zeros(size(PSF));
while S(i,i)/S(1,1) >= Rtol
  a = R \ ( U(:,i) * sqrt(S(i,i)) );
  b = R \ ( V(:,i) * sqrt(S(i,i)) );
  PP = PP + a*b';
  % Comment out the next two statements, replace with the two that
  % follow, in order to fix the 11/22/02 problem.
  %A{i} = build_toep(a, center(1), m) + buildHank(a, center(1), m);
  %B{i} = build_toep(b, center(2), n) + buildHank(b, center(2), n);
  B{i} = build_toep(a, center(1), m) + buildHank(a, center(1), m);
  A{i} = build_toep(b, center(2), n) + buildHank(b, center(2), n);
  i=i+1;
end


%
%  The following scales the approximation so that the
%  corresponding PSF approximation has the property
%  that its sum of values is the same as the sum of the
%  original PSF values.
%
cs = sqrt(sum(PSF(:))/sum(PP(:)));
for i = 1:length(A)
  A{i} = cs*A{i};
  B{i} = cs*B{i};
  AA = A{i};
  BB = B{i};
end


