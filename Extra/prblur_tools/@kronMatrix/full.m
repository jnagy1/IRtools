function A = full(K)
%
%                A = full(K);
%
%  Converts a kronMatrix matrix K to full storage organization.  
%
A = kron(K.a{1}, K.b{1});
for i = 2:length(K)
  A = A + kron(K.a{i}, K.b{i});
end
  