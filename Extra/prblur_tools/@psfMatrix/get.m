function val = get(A, prop_name)
%
%  Get psfMatrix properties from specified object and return
%  the value.
%
%     val = get(A, prop_name);
%
%  Valid choices for prop_name are:
%    'psf', 'matdata', 'type', 'boundary', 'transpose'
%

%  J. Nagy  5/2/01

switch prop_name
case 'psf'
  val = A.psf;
case 'matdata'
  val = A.matdata;
case 'type'
  val = A.type;
case 'boundary'
  val = A.boundary;
case 'transpose'
  val = A.transpose;
otherwise
  error([prop_name, 'Is not valid psfMatrix property'])
end