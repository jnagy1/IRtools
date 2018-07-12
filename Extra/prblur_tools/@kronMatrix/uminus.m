function m = uminus(k)
%  kronMatrix/uminus
%
%     -k negates the elements of k.
%
%     -kron(k) == kron(-k)

m = k;
if (length(m.a) > length(m.b))
   m.b = -(m.b);
else
   m.a = -(m.a);
end

