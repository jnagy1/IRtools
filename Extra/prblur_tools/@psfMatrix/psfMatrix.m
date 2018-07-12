function A = psfMatrix(varargin)
%
%  CONSTRUCTOR FOR psfMatrix OBJECT
%
%  This new version uses Kronecker products, and should correct a bug
%  in the old version, which gave incorrect results for A'*x with
%  a non-symmetric PSF and refexive boundary conditions.
%
%  The psfMatrix class is based on a structure with six fields:
%    psf      - psf object 
%    matdata  - matrix data needed to do matrix vector multiplications
%    transp_matdata - needed for multiplying by transpose.
%    type     - character string indicating:
%                 'invariant', 'variant', 'separable'
%    boundary - character string array indicating type of boundary conditions
%               to be used.  Choices are:
%                 'zero', 'periodic', 'reflexive' (or 'neumann')
%               The default is reflexive.
%    transpose- indicates if the matrix has been transposed.
%    imsize   - indicates the size of the image.  This might be
%               needed for some space variant problems.
%
%  Calling Syntax:
%       A = psfMatrix             (returns object with empty fields)
%       A = psfMatrix(psfMatrixObj) 
%       A = psfMatrix(PSF, boundary, center, imsize)
%
%    where 
%       * psfMatrixObj is an already existing psfMatrix object
%       * psfObj       is a psf object (see help psf for more information)
%       * PSF          can be either a double array containing one PSF
%                      image, or a cell array containing one or more 
%                      PSF images
%       * boundary     is a character string indicating desired boundary
%                      condition. (see above)
%       * center       is either a double array or a cell array containing
%                      the (i,j) locations of the point sources of the 
%                      PSF images
%       * imsize       size of reconstructed image in deblurring problem
%                      While it is not necessary to give this value, it 
%                      can improve speed of iterative methods because it
%                      allows pre-computation of certain terms
%

%  J. Nagy & K. Lee  1/12/02
%  11/20/06 -- made reflexive default BC
%  7/5/16 fixed reflexive boundary conditions transpose mult
%

switch nargin

case 0
  A.psf = psf;
  A.matdata = [];
  A.transp_matdata = [];
  A.type = '';
  A.boundary = '';
  A.transpose = 0;
  A.imsize = [];
  A = class(A, 'psfMatrix');
  return

case 1
  if ( isa( varargin{1}, 'psfMatrix' ) )
    A = varargin{1};
    return
  end
  PSF = varargin{1};
  boundary = [];
  center = [];
  imsize = [];
case 2
  PSF = varargin{1};
  boundary = varargin{2};
  center = [];
  imsize = [];
case 3
  PSF = varargin{1};
  boundary = varargin{2};
  center = varargin{3};
  imsize = [];
case 4
  PSF = varargin{1};
  boundary = varargin{2};
  center = varargin{3};
  imsize = varargin{4};
end
if isempty(boundary), boundary = 'reflexive'; end
if isempty(imsize)
    if isa(PSF, 'double')
        imsize = size(PSF);
    end
end
if isa(PSF, 'cell')
    A.type = 'variant';
else
    A.type = 'invariant';
end
if isempty(center)
    P = psf(PSF);
else
    P = psf(PSF, center);
end
A.psf = P;
A.boundary = boundary;
if strcmp(A.type, 'variant')
    A.matdata = constructMatrix( P.image, P.center );
    A.transp_matdata = constructTranspMatrix(P.image, P.center);
else
    centerP = P.center;
    center = centerP{1};
    A.center = center;
    A.psf = P;
    [A.matdata, A.p] = constructMatrix_new( PSF, center, boundary, imsize );
end
A.transpose = 0;
A.imsize = imsize;
A = class(A, 'psfMatrix');
