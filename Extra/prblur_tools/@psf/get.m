function val = get(P, prop_name)
%
%  GET Get psf properties from specified object and return
%      the value.
%
%     val = get(P, prop_name);
%
%  Valid choices for prop_name are:
%    'image' and 'center'
%

%  J. Nagy 5/1/01

switch prop_name
case 'image'
  val = P.image;
case 'center'
  val = P.center;
otherwise
  error([prop_name, 'Is not valid psf property'])
end