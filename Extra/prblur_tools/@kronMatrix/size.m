function varargout = size(k)
%  kronMatrix/size
%
%     returns the size of a kronMatrix object.
%

s = size(k.a) .* size(k.b);
switch (nargout)
  case 0
    disp(sprintf('\nans =\n'));
    disp(s);
  case 1
    varargout{1} = s;
  case 2
    varargout{1} = s(1);
    varargout{2} = s(2);
end
