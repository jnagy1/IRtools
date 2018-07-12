function n = kron_mtimes(k, m)

%  helper function for kronMatrix mtimes
%
%
eval('n = kronMatrix(k.a * m.a, k.b * m.b);', ...
     'error(lasterr)')
end
