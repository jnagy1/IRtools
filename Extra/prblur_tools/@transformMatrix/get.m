function val = get(A, prop_name)
%
%  Get transformMatrix properties from specified object and return
%  the value.
%
%     val = get(A, prop_name);
%
%  Valid choices for prop_name are:
%    'transform', 'transpose'
%

%  J. Nagy  6/2/02

switch prop_name
case 'transform'
  val = A.transform;
case 'transpose'
  val = A.transpose;
otherwise
  error([prop_name, 'Is not valid transformMatrix property'])
end