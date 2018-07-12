function G = kronMatrix( varargin )

%
% G = kronMatrix( varargin )
%
% Constructor for the kronMatrix class.  This class is used to
% represent objects of the form:
%       \sum [A{i} (x) B{i}]
%
% Object's fields are:
%    G.a = cell array containing A{i} if i>1.
%    G.b = cell array containing B{i} if i>1.          
%
% Calling Syntax:
%       G = kronMatrix
%       G = kronMatrix(G);
%       G = kronMatrix(A, B);
%
%   where
%       G = a pre-existing kronMatrix object
%       A, B = both matrices or both cell arrays (same length) such that
%              \sum [A{i} (x) B{i}] represents a bigger matrix.
%

%  L. Perrone, 5/13/03

%  Modifications:
%  6/4/03, J. Nagy
%          Cosmetic changes to incorporate into RestoreTools
%

%
% build the kronMatrix based on number and type of input arguments
%
switch nargin

  case 0
   G.a = cell(0);
   G.b = cell(0);
   G = class(G, 'kronMatrix');
 
  case 1
   if isa(varargin{1}, 'kronMatrix')
      G = varargin{1};
   else
      error('Wrong input argument')
   end

  case 2

   if isa(varargin{1}, 'double') & isa(varargin{2}, 'double')
     A = cell(1);, A{1} = varargin{1};
     B = cell(1);, B{1} = varargin{2};
     G.a = A;
     G.b = B;
     G = class(G, 'kronMatrix');

   elseif isa(varargin{1}, 'cell') & isa(varargin{2}, 'cell')
     if length(varargin{1}) == length(varargin{2})     
        G.a = varargin{1};
        G.b = varargin{2};
        G = class(G, 'kronMatrix');
     else
        error('Input cell arrays must have the same length')
     end

   else
     error('Wrong input arguments')

   end % case 2

  otherwise
    error('Wrong input arguments')

end % switch nargin