function B = subsref(A, index)
%
%  Define field name indexing for transformMatrix object.
%
%               B = subsref(A, index)
%
%          This is called whenever a subscripted reference to the
%          transformMatrix object is made, such as:
%              A.fieldname, fieldname = transform, transpose
%

%  J. Nagy  1/18/02

switch index.type
case '()'
  error('Paranthetical indexing not supported for transformMatrix object.')

case '.'
  switch index.subs
  case 'transform'
    B = A.transform;
  case 'transpose'
    B = A.transpose;
  otherwise
    error('Invalid field name.')
  end

case '{}'
  error('Cell array indexing not supported for transformMatrix object.')
end
