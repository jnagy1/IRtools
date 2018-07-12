function c=vertcat(varargin)
%
%           vertical concatination of the psfMatrix objects
%
%  Input:  any number of psfMatrix objects, or a multiPsfMatrix object 
%
%  Output:  c is a multiPsfMatrix object of size (number of argin) by 1
%

%  K. Lee  2/28/02

if isa( varargin{1}, 'cell')
  A=varargin{1};  			   % A is a cell
  A=A{1};         			   % A is a psfMatrix object
  bdry=A.boundary;
  Apsf=A.psf;        			   % Apsf is a cell
  [r,d]=imageSize(A);
else
  A=varargin{1};
  bdry=A.boundary;
  Apsf=A.psf;
  [r,d]=imageSize(A);
end
c=multiPsfMatrix(varargin{1});
for i=2:nargin
  A=varargin{i};
  [r1,c1]=imageSize(A);
  bdry2=A.boundary;
    if ~isa(varargin{i},'psfMatrix'),
      error('psfMatrix objects can only be concatenated with other psfMatrix objects')
  elseif  (strcmp(bdry,bdry2)==1) 
        c=multiPsfMatrix(c,varargin{i});
    else
      error('psfMatrix objects must have the same boundary')

    end
end