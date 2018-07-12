function [a,b,c,d] = parsesubs(x, subs)
%
%


if (subs{1} == ':')
   s = size(x);
   a = 1;
   b = s(1);
elseif (isa(subs{1}, 'double'))
   a = subs{1};
   b = subs{1};
end

c = 1;
d = 1;