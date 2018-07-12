function b = mtimes(A, x)
%  
%   Overload matrix multiplication operation for transformMatrix
%
%   Note:  We scale fft so that the transformMatrix is unitary.
%     

%  J. Nagy, 6/2/02

n = prod(size(x));

if (isa(A, 'transformMatrix'))
  switch A.transform
  case 'fft'
    if A.transpose
      b = sqrt(n) * ifftn(x);
    else
      b = fftn(x) / sqrt(n);
    end
    if ( max(abs(imag(b(:)))) <= sqrt(eps) )
      b = real(b);
    end
  case 'dct'
    if A.transpose
      if (ndims(x) == 1)
        b = idct(x);
      elseif (ndims(x) == 2)
        b = idct2(x);
      else
        error('Can only do this for 1 and 2 dimensional vectors')
      end
    else
      if (ndims(x) == 1)
        b = dct(x);
      elseif (ndims(x) == 2)
        b = dct2(x);
      else
        error('Can only do this for 1 and 2 dimensional vectors')
      end
    end
  otherwise
    error('Incorrect transform type.')
  end
else
  error ('Argument 1 must be a transformMatrix')
end

