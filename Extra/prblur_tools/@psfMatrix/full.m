function Afull = full( A )
%
%  Overload full for psfMatrix object.
%  Be carefull not to use this for large images -- maybe
%  max size should be 64-by-64.
%

% J. Nagy  1/23/2016

K = kronApprox(A, [], 0);
Afull = zeros(size(A));
if iscell(K.a)
    kron_terms = length(K.a);
else
    kron_terms = 1;
end
for j = 1:kron_terms
    Afull = Afull + kron(K.a{j}, K.b{j});
end