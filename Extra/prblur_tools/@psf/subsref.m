function B = subsref(P, index)
%
%  SUBSREF Define field name indexing for psf object.
%
%               B = subsref(P, index)
%
%          This is called whenever a subscripted reference to the
%          psf object is made, such as:
%              P(i), i = 1, 2
%              P.fieldname, fieldname = image, center
%

%  J. Nagy  5/1/01

switch index.type
case '()'
  switch index.subs{:}
  case 1
    B = P.image;
  case 2
    B = P.center;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'image'
    B = P.image;
  case 'center'
    B = P.center;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for psf object.')
end
