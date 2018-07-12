function l = length(K)
%
%  Returns the number of terms in the sum:
%       K = \sum [ A{i} (x) B{i} ]
%
%  where K is a kronMatrix object.
%
if isa(K.a, 'cell')
  l = length(K.a);
else
  l = 1;
end