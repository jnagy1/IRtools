function A = subsasgn(A, index, val)
%
%  SUBSASGN  Define index assignment for transfromMatrix object.
%
%               A = subsref(A, index, val)
%
%          This is called whenever an assignment statement of a
%          transformMatrix object is made, such as:
%              A.fieldname = val, fieldname = transform, transpose
%

%  J. Nagy 6/2/02

switch index.type
case '()'
  error('Parenthetical indexing not supported for transformMatrix object')

case '.'
  switch index.subs
  case 'transform'
    P.transform = val;
  case 'transpose'
    P.transpose = val;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for transformMatrix object.')
end
