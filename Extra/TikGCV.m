function G = TikGCV(alpha, bhat, s, nORomega)
%
%    G = TikGCV(alpha, bhat, s, nORomega)
%
%  This function evaluates the GCV function for Tikhonov
%  regularization in standard form.  
%
%  Input:  alpha -  regularization parameter
%           bhat -  vector U'*b, where U = left singular vectors
%                   or left generalized singular vectors
%              s -  vector containing the singular values or
%                   generalized singular values
%       nORomega -  used to implement either the weighted GCV function
%                   (with 0<omega<=1) or the GCV function with explicit 
%                   appearence of the 'degrees of freedom')
%
%  Output:     G -  the scalar G(alpha).
%
%  J.Chung and J. Nagy  3/2007
%  Updated by S.Gazzola 7/2017

if nargin == 3
  nORomega = [];
end
if isempty(nORomega)
  nORomega = 1;
end

m1 = length(bhat);
m2 = length(s);

if nORomega <= 1
    wGCV = 1;
    omega = nORomega;
else
    wGCV = 0;
    n = nORomega;
    dof = n - m2;
end

t0 = sum(abs(bhat(m2+1:m1)).^2);

s2 = abs(s) .^ 2;
alpha2 = alpha^2;

tt = 1 ./ (s2 + alpha2);
t1 = alpha2 .* tt;
t2 = abs(bhat(1:m2) .* t1) .^2;

if wGCV
    t3 = t1 + (1 - omega)*s2 .* tt;
    G = (sum(t2) + t0) / (m1 - m2 + sum(t3))^2;
else
    G = (sum(t2) + t0) / (dof + sum(t1))^2;
end
