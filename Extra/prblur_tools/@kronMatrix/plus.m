function n = plus(k,m)
%  kronMatrix +
%
%     returns the sum of two kronMatrix objects.
%
%     n = k + m will result in a kronMatrix object if
%        k.a == m.a or k.b == m.b
%
%     otherwise, this function will return an explicitly
%        formed matrix, i.e. kron(k) + kron(m)
%

if (isa(k, 'kronMatrix'))
   if (isa(m, 'kronMatrix'))
      if isequal(k.a, m.a)
         eval('n = kronMatrix(k.a, k.b + m.b);', ...
              'error(lasterr)')
      elseif isequal(k.b, m.b)
         eval('n = kronMatrix(k.a + m.a, k.b);', ...
              'error(lasterr)')
      else
         eval ('n = kron(k) + kron(m);', ...
               'error(lasterr)')
      end
   else
      eval('n = kron(k) + m;', ...
           'error(lasterr)')
   end
else
   eval('n = kron(m) + k;', ...
        'error(lasterr)')
end

