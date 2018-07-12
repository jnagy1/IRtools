function P = subsasgn(P, index, val)
%
%  SUBSASGN  Define index assignment for psf object.
%
%               P = subsref(P, index, val)
%
%          This is called whenever an assignment statement of a
%          psf object is made, such as:
%              P(i) = val, i = 1, 2
%              P.fieldname = val, fieldname = image, center
%

%  J. Nagy 5/1/01

switch index.type
case '()'
  switch index.subs{:}
  case 1
    P.image = val;
  case 2
    P.center = val;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'image'
    P.image = val;
  case 'center'
    P.center = val;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for psf object.')
end
