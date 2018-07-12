function B = subsref(A, index)
%
%  Define field name indexing for psfMatrix object.
%
%               B = subsref(A, index)
%
%          This is called whenever a subscripted reference to the
%          psfMatrix object is made, such as:
%              A(i)
%              A.fieldname, fieldname = psf, matdata, transp_matdata, type, boundary, transpose
%

%  J. Nagy  1/18/02

switch index.type
case '()'
  switch index.subs{:}
  case 1
    B = A.psf;
  case 2
    B = A.matdata;
  case 3
    B = A.transp_matdata;
  case 4
    B = A.type;
  case 5
    B = A.boundary;
  case 6
    B = A.transpose;
  case 7
    B = A.imsize;
  otherwise
    error('Index out of range.')
  end
  
case '.'
  switch index.subs
  case 'psf'
    B = A.psf;
  case 'matdata'
    B = A.matdata;
  case 'transp_matdata';
    B = A.transp_matdata;
  case 'type'
    B = A.type;
  case 'boundary'
    B = A.boundary;
  case 'transpose'
    B = A.transpose;
  case 'imsize'
    B = A.imsize;
  case 'p'
    B = A.p;
  case 'center'
    B = A.center;
  otherwise
    B = subsref(A.psf, index);
  end

case '{}'
  error('Cell array indexing not supported for psfMatrix object.')
end
