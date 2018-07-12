function nrm = norm( A, p )
%
%  Overload norm for psfMatrix object.
%  This will only compute a 1 or infinity norm.
%
%  Note that here we assume that A has no negative values,
%  so the 1 or infinity norm can be computed by matrix-vector
%  products:
%     norm(A,inf) = max(A*ones(n,1))
%     norm(A,1) = max(A'*ones(m,1))
%

% J. Nagy  8/17/11

if nargin == 1
  p = 1;
  warning('Computing default 1-norm')
end
[m, n] = size(A);
switch p
  case inf
    nrm = max(A*ones(n,1));
  case 'inf'
    nrm = max(A*ones(n,1));
  case 1
    nrm = max(A'*ones(m,1));
  otherwise
    error('Can only compute 1 or infinity norm of a psfMatrix')
end
