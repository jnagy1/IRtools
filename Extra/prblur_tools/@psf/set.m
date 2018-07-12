function P = set(P, varargin)
%
% SET Set psf properties and return the updated object.
%
%     P = set(P, ...)
% 
%  This funtion accepts a psf object, P, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'image', 'center'
%

%  J. Nagy 5/2/01

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'image'
    P.image = val;
  case 'center'
    P.center = val;
  otherwise
    error('Valid psf properties: image, center')
  end
end
