function out = mpower(k, in2)

%  kronMatrix k^x
%
%  raises a kronMatrix object to a power.
%
%  k = a (x) b ==> k^x = (a (x) b)^x

if (isa(k, 'kronMatrix') & isa(in2, 'double') & size(in2) == 1 & prod(size(k)) == length(k)^2)
   out = kronMatrix(k.a^in2, k.b^in2);
else error('k must be a square kronMatrix and x must be a scalar');
end

