
function [A, B] = zeroKronApprox( K, m, n, Rtol)
%
%       [A, B] = zeroKronApprox( K, m, n, Rtol);
%
% Helper function for kronApprox.m;
% computes a kronecker sum approximation to the blurring matrix K
% that arises from the input PSF under zero boundary conditions.
% This approximation is done a-la Kamm-Nagy (see paper for details).  
%
%  Input:
%         K - psfMatrix object
%         m - size of matrix A (assumed square)
%         n - size of matrix B (assumed square)
%         Rtol - relative tolerance ( < 1 ); determines # of kron products
%                in the approximation.  If si/s1 > Rtol, 
%                where si is the ith largest
%                singular value of the weighted PSF,
%                then Ai(x)Bi will be included in the approximation.
%                OR
%                positive integer specifying the number of terms in the
%                approximation.
%
%  Output:
%      Cells A and B such that K \approx  \sum( A{i} (x) B{i} )
%

%  J. Nagy  2/11/02

%  Modifications:  This used to be zeroSepMatrix, which did not compute
%                  A and B. 
%  5/25/02, J. Nagy 
%           This now computes A and B.  
%  9/02     L. Perrone
%           This code now conforms to the new kronMatrix class,
%           where fields K.a and K.b are cell arrays containing
%           (possibly) more than one matrix.
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
A = cell(1);
B = cell(1);
PP = zeros(size(PSF));
if Rtol >= 1 && Rtol == fix(Rtol)
    MaxTerms = Rtol;
    MaxTerms0 = MaxTerms;
    Rtol = 0;
else
    MaxTerms = min(min(size(S)));
    MaxTerms0 = -1;
end
Rtol = Rtol + 10*eps;
MaxTerms = min(sum(diag(S)/S(1,1) >= Rtol), MaxTerms);
if MaxTerms0 > 0 && MaxTerms0 ~= MaxTerms
    warning('Max number of terms in Kronecker sum is only %d', MaxTerms)
end

for i = 1:MaxTerms
  a = ( U(:,i) * sqrt(S(i,i)) ) ./ Aweights;
  b = ( V(:,i) * sqrt(S(i,i)) ) ./ Bweights;
  PP = PP + a*b';
  % Comment out the next two statements b/c they rely on vec ordering
  % of the PSF, not lexicographical... 
  % Replace with the two statements that follow:
  %  A{i} = build_toep(a, center(1), m);
  %  B{i} = build_toep(b, center(2), n);
  B{i} = build_toep(a, center(1), m);
  A{i} = build_toep(b, center(2), n);
  % now resume business as usual
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