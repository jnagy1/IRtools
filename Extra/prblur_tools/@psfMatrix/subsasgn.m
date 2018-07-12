function P = subsasgn(P, index, val)
%
%  SUBSASGN  Define index assignment for psfMatrix object.
%
%               P = subsref(P, index, val)
%
%          This is called whenever an assignment statement of a
%          psf object is made, such as:
%              P(i) = val, i = 1, 2, 3, 4, 5, 6
%              P.fieldname = val, 
%                fieldname = psf, matdata, transp_matdata, type, boundary, transpose, imsize
%

%  J. Nagy 5/1/01

switch index.type
case '()'
  switch index.subs{:}
  case 1
    P.psf = val;
  case 2
    P.matdata = val;
  case 3
    P.transp_matdata = val;
  case 4
    P.type = val;
  case 5
    P.boundary = val;
  case 6
    P.transpose = val;
  case 7
    P.imsize = val;
  otherwise
    error('Index out of range.')
  end

case '.'
  switch index.subs
  case 'psf'
    P.psf = val;
  case 'matdata'
    P.matdata = val;
  case 'transp_matdata'
    P.transp_matdata = val;
  case 'type'
    P.type = val;
  case 'boundary'
    P.boundary = val;
  case 'transpose'
    P.transpose = val;
  case 'imsize'
    P.imsize = val;
  otherwise
    error('Invalid field names.');
  end

case '{}'
  error('Cell array indexing not supported for psf object.')
end
