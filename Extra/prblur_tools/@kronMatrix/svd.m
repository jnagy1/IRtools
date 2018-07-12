function varargout = svd(G)
%
% [U,S,V] = svd(G);
% s = svd(G);
%
% Computes an approximate (usually) svd of the matrix represented
% by the kronMatrix object G = \sum [ A{i} (x) B{i} ].  If there
% is only one term in the sum, then the svd is "exact".
%
% On entry: G = a kronMatrix object
% On exit: U, V, S = kronMatrix objects 
%          S = vector of singular values (not in descending order b/c
%              of the kron business)
%            ... such that U*S*V' is an approximation to G

%
% First get U and V
%
A1 = G.a{1};
B1 = G.b{1};
[Ua, Sa1, Va] = svd(A1);
[Ub, Sb1, Vb] = svd(B1);

%
% Now work through the computations to get S (Kamm-Nagy LAA)
%
l = length(G);
Sa = cell(l, 1);
Sb = cell(l, 1);
Sa{1} = Sa1;
Sb{1} = Sb1;
for i = 2:l
  Sa{i} = diag(diag(Ua' * G.a{i} * Va));
  Sb{i} = diag(diag(Ub' * G.b{i} * Vb));
end

if (nargout == 3)
  varargout{1} = kronMatrix(Ua, Ub);
  varargout{2} = kronMatrix(Sa, Sb);
  varargout{3} = kronMatrix(Va, Vb);
elseif (nargout == 1)
  s = kron(diag(Sa{1}), diag(Sb{1}));
  for i = 2:length(G);
    s = s + kron(diag(Sa{i}), diag(Sb{i}));
  end
  s = flipud(sort(s));
  varargout{1} = s;
else
  error('Incorrect number of output arguments.')
end

