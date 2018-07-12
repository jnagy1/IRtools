function varargout = svdKron( A, b, k, kronTol )
%
%       [U, s, V] = svd(A);
%       [U, s, V] = svd(A, b);
%       [U, s, V] = svd(A, b, k);
%       [U, s, V] = svd(A, b, k, kronMax);
%
%  Overloaded SVD method for psfMatrix objects. This computes the SVD
%  componets for a fixed rank.
%
%  Input: 
%    A is a psfMatrix
%
%  Optional Input:
%    b is a blurred image.  Since A often uses a compact storage
%      scheme, this is sometimes needed to determine "real size" 
%      of the matrix.
%    k is a positive integer indicating the number of largest singular 
%      values and corresponding vectors to compute. 
%      Default value is k = 10.
%    kronTol   relative tolerance ( < 1 ); determines number of terms in
%              in the sum of Kronecker product decomposition.  Specifically,
%              If si/s1 > Rtol, where si is the ith largest
%              singular value of the weighted PSF,
%              then Ai(x)Bi will be included in the approximation.
%                OR
%              positive integer specifying the maximum number of terms in the
%              approximation.
%              Default value is kronTol = 1.
%
%  Output:
%

%  J. Nagy  1/21/16


%switch A.type
%  case 'invariant'
%    P = A.psf;
%    P1 = P.image;, c = P.center;
%    PSF = P1{1};,  center = c{1};
%
%    if (nargin == 2)
%      PSF = padarray(PSF, size(b) - size(PSF), 'post');
%    else
%      b = PSF;  % This is used to define correct dimensions only.
%    end
%  
%    switch A.boundary
%    case 'periodic'
%      [U, S, V] = fft_svd(PSF, center);
%    case 'reflexive'
%      if issymmetric(A)
%         [U, S, V] = dct_svd(PSF, center);
%      else
%         [U, S, V] = svd_approx( kronApprox(A, b) );
%      end
%    case 'zero'
%      [U, S, V] = svd_approx( kronApprox(A, b) );
%    otherwise
%      error('Invalid boundary condition')
%    end
%   
%  case 'variant'
%    [U, S, V] = svd_approx( kronApprox(A, b) );
%    
%  otherwise
%    error('invalid matrix type')
%end
if nargin == 1
    b = []; k = []; kronTol = [];
elseif nargin == 2
    k = []; kronTol = [];
elseif nargin == 3
    kronTol = [];
end
if isempty(k), k = 10; end
if isempty(kronTol), kronTol = 1; end

K = kronApprox(A, b, kronTol);

%
%  Get first term of the Kronecker sum, and compute its full SVD.
%
if iscell(K.a)
    Ka1 = K.a{1};
    Kb1 = K.b{1};
    kronTerms = length(K.a);
else
    Ka1 = K.a;
    Kb1 = K.b;
    kronTerms = 1;
end
%K1 = kronMatrix(Ka1,Kb1);
%[U1, S1, V1] = svd(K1);
[U1a, S1a, V1a] = svd(Ka1);
[U1b, S1b, V1b] = svd(Kb1);

%
%  Now compute the reduced SVD, based on the given rank, k.
%
%U1a = U1.a; U1b = U1.b;
%V1a = V1.a; V1b = V1.b;
%
%  This is a stupid hack -- need to fix this later.  I want to grab
%  at least the largest k values, but I have to grab some from each
%  Kronecker term.  So I may return more than the specified k.
%
%s_all = kron(diag(S1.a{1}),diag(S1.b{1}));
s_all = kron(diag(S1a), diag(S1b));
[~, idx_all] = sort(s_all, 'descend');
idx = idx_all(1:k);
n = min(size(Ka1));
Aidx = ceil(idx/n); 
Bidx = idx - (Aidx-1).*n;
Aidx_max = max(Aidx);
Bidx_max = max(Bidx);
U1 = kronMatrix(U1a, U1b);
V1 = kronMatrix(V1a, V1b);
U1a = U1a(:,1:Aidx_max); U1b = U1b(:,1:Bidx_max);
V1a = V1a(:,1:Aidx_max); V1b = V1b(:,1:Bidx_max);
k_new = Aidx_max*Bidx_max;
if k ~= k_new
    warning('could not give specified rank %d, increasing to rank %d',k,k_new)
end
T = zeros(k_new, k_new);
for j = 1:kronTerms
  T = T + kron(U1a'*K.a{j}*V1a,U1b'*K.b{j}*V1b);
end
[Ut, S, Vt] = svd(T);
% [min(diag(S)); max(s_all(k_new+1:end))]
s = [diag(S); s_all(k_new+1:end)];

if nargout == 1
  varargout{1} = 2;
elseif nargout == 3
  U.kron = U1; U.small = Ut;
  V.kron = V1; V.small = Vt;
  varargout{1} = U;
  varargout{2} = s;
  varargout{3} = V;
else
  error('In correct number of output variables')
end