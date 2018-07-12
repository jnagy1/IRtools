function N = mtimes(K, M)
% N = mtimes(K, M)
%
%  kronMatrix multiplication;
%     multiply a kronMatrix by a matrix, a vector, or by
%     another kronMatrix (if possible),
%     
%     

% 9/2002 L. Perrone 
% written for new kronMatrix class

% 6/2003 J. Nagy
% some minor modifications to incorporate new subsref.m capabilities

if (isa(K, 'kronMatrix'))
   if isa(M, 'double')
      N = left_mtimes(K, M);
   elseif isa(M,'kronMatrix')
      if length(M) == 1
         sizeMa = size(M.a{1});
         sizeMb = size(M.b{1});
         l = length(K);
         Anew = cell(l,1);
         Bnew = cell(l,1);
         for i = 1:l
           sizeKa = size(K.a{i});
           sizeKb = size(K.b{i});
           if sizeKa(2)==sizeMa(1) & sizeKb(2)==sizeMb(1)
             Anew{i} = K.a{i} * M.a{1};
             Bnew{i} = K.b{i} * M.b{1};
           else
             error('Kron factors must be of compatible sizes for multiplication')
           end
         end % end the for i=1:l loop
         N = kronMatrix(Anew,Bnew);
      elseif length(K) == 1
         sizeKa = size(K.a{1});
         sizeKb = size(K.b{1});
         l = length(M);
         Anew = cell(l,1);
         Bnew = cell(l,1);
         for i = 1:l
           sizeMa = size(M.a{i});
           sizeMb = size(M.b{i});
           if sizeKa(2)==sizeMa(1) & sizeKb(2)==sizeMb(1)
             Anew{i} = K.a{1} * M.a{i};
             Bnew{i} = K.b{1} * M.b{i};
           else
             error('Kron factors must be of compatible sizes for multiplication')
           end
         end % end the for i=1:l loop
         N = kronMatrix(Anew,Bnew);
      else 
         error('This currently only works if one of the kronMatrix objects has one term in it''s sum.')
      end
   else % M isn't a 'double' or a 'kronMatrix'
      error('Wrong input arguments')
   end % end the case when K is a 'kronMatrix'

elseif isa(K,'double') & isa(M,'kronMatrix')
   N = right_mtimes(K, M);
else
   error('Wrong input arguments')
end

