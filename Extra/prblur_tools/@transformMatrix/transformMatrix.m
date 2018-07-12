function A = transformMatrix (varargin)
%  
%  CONSTRUCTOR FOR transfromMatrix OBJECT
%
%  The transformMatrix class is based on a structure with the
%  two fields:
%    transform  - character string indicating the transform:
%                 'fft'  or 'dct'
%    transpose  - indicates if the matrix has been transposed.
%
%   Calling Syntax:
%       A = transformMatrix           (returns object with empty fields)
%       A = transformMatrix(transformMatrixObj)
%       A = transformMatrix(transform)
%
%    where
%       * transformMatrixObj is an already existing transformMatrix object
%       * transform is a character string indicating the transform:
%         'fft' or 'dct'
%

%  J. Nagy  6/2/02


switch nargin
case 0
% default
   A.transform = '';
   A.transpose = 0;
   A = class(A, 'transformMatrix');
case 1
% 
   if isa (varargin{1}, 'kronMatrix')
      A = varargin{1};
   elseif ischar(varargin{1})
      A.transform = varargin{1};
      A.transpose = 0;
      A = class(A, 'transformMatrix');
   end
otherwise
   error ('Incorrect number of arguments');
end

