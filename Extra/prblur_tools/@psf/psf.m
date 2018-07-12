function P = psf(varargin)
%
%  CONSTRUCTOR FUNCTION FOR PSF OBJECT
%
%  The psf class is based on a structure with two fields:
%    image - cell array containing the images of one or more
%            point spread functions.
%
%    center - cell array containing the centers of the individual
%             PSFs (that is, the locations of the point sources
%             which created the PSFs).
%
%  Calling Syntax:
%       P = psf
%       P = psf(psfObj)
%       P = psf(PSF)
%       P = psf(PSF, center)
%
%  where
%    * psfObj   is already a psf object
%    * PSF      can either be a double array containing a single PSF
%               image, or a cell array containing one or more PSF images.
%    * center   is either a double array of length two, or a cell array
%               containing the center(s) of the PSF(s).
%               If the center is not specified, we assume it's at or
%               near the location of the maximum entry of the PSF.
%

%  J. Nagy 1/12/02

switch nargin
case 0
  % if no input arguments, create default object
  P.image = {};
  P.center = {};
  P = class(P, 'psf');
case 1
  % if single argument of class psf, return it
  if ( isa( varargin{1}, 'psf' ) )
    P = varargin{1};
  elseif ( isa (varargin{1}, 'double' ) )
    I = cell(1);
    I{1} = varargin{1};
    P.image = I;
    P.center = findCenter( P.image );
    P = class(P, 'psf');
  elseif ( isa (varargin{1}, 'cell' ) )
    P.image = varargin{1};
    P.center = findCenter( P.image );
    P = class(P, 'psf');
  else
    error('Incorrect argument type')
  end

case 2
  % create object using specific values
  if ( isa( varargin{1}, 'double' ) )
    I = cell(1);
    I{1} = varargin{1};
    P.image = I;
    if ( isa( varargin{2}, 'double' ) )
      c = cell(1);
      c{1} = varargin{2};
      P.center = c;
    elseif ( isa( varargin{2}, 'cell') )
      P.center = varargin{2};
    else
      error('Incorrect input arguments')
    end
  elseif ( isa( varargin{1}, 'cell' ) )
    P.image = varargin{1};
    if ( isa( varargin{2}, 'cell' ) )
      if ( prod(size(varargin{1})) == prod(size(varargin{2})) )
         P.center = varargin{2};
      else
         error('Input arguments incorrect')
      end
    else
      error('Input arguments incorrect')
    end
  end
  if isempty(varargin{2})
      P.center = findCenter(P.image);
  end
  P = class(P, 'psf');
otherwise
  error('Incorrect number of input arguments')
end

