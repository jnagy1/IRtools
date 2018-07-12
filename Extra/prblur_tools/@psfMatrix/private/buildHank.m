function H = buildHank(c, k, n)

%     H = buildHank(c, k, n);
%
% builds a hankel matrix of dimension n, around b, k as in lpjn paper
%

% 4/26/02, L. Perrone.
%
% Modifications:
% 5/25/01, J. Nagy
%          Some cosmetic changes to incorporate into RestoreTools
%

m = length(c);

col = zeros(n,1);
col(1:m-k) = c(k+1:m);

row = zeros(n,1);
row(n-k+2:n) = c(1:k-1);

H = hankel(col, row);