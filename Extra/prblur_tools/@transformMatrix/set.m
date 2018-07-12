function A = set(A, varargin)
%
% Set transformMatrix properties and return the updated object.
%
%     A = set(A, ...)
% 
%  This funtion accepts an transformMatrix object, A, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'transform', 'transpose'
%

%  J. Nagy  6/2/02

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'transform'
    A.transform = val;
  case 'transpose'
    A.transpose = val;
  otherwise
    error('Valid transformMatrix properties: transform, transpose')
  end
end
  