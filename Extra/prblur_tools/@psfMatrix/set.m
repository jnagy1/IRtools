function A = set(A, varargin)
%
% Set psfMatrix properties and return the updated object.
%
%     A = set(A, ...)
% 
%  This funtion accepts an psfMatrix object, A, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'psf', 'matdata', 'type', 'boundary', 'transpose'
%

%  J. Nagy  5/2/01

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'psf'
    A.psf = val;
  case 'matdata'
    A.matdata = val;
  case 'type'
    A.type = val;
  case 'boundary'
    A.boundary = val;
  case 'transpose'
    A.transpose = val;
  otherwise
    error('Valid psfMatrix properties: psf, matdata, type, boundary, transpose')
  end
end
  