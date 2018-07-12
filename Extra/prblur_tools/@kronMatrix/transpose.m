function out = transpose(in)

%  kronMatrix k.'
%    returns the transpose of the kronMatrix object k
%
%    k = \sum [a{i} (x) b{i}] ==> k.' = \sum [a{i}.' (x) b{i}.']
%
%

% Spring 2002 created by R. Wright
%

% 9/2002 L. Perrone 
% Modified code to bring it into line with the new kronMatrix class
%

A = in.a;
B = in.b;
l=length(A);
for i=1:l
 tmp = A{i};
 A{i} = tmp.';
 tmp = B{i};
 B{i} = tmp.';
end
out = kronMatrix(A,B);