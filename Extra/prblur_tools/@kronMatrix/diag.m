function m = diag(k)
%  kronMatrix diag(k)
%
%     returns a kronMatrix object m, such that
%     m.a contains the diagonal elements of k.a
%     and m.b contains the diagonal elements of k.b
%
%     kron(diag(k)) = diag(kron(k)).
%

m = kronMatrix(diag(k.a),diag(k.b));

