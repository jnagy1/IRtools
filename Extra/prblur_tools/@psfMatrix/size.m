function varargout = size( A, dim )
%
%  Overload size for psfMatrix object.
%  This just displays the storage requirements for the matrix
%  data.
%

% J. Nagy  1/9/02

if nargin == 2
  d = size(A.psf);
  if dim <= length(d)
    d = d(dim);
  else 
    d = 1;
  end
else
  d = size(A.psf);
end
d = d .^2;


if nargout == 1 || nargout == 0
  varargout{1} = d;
else
  for i = 1:length(d)
    varargout{i} = d(i);
  end
end