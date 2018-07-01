function [A, b, x, ProbInfo] = PRspherical(varargin) 
% PRspherical Generates data for spherical means tomography problems
%
% [A, b, x, ProbInfo] = PRspherical
% [A, b, x, ProbInfo] = PRspherical(n)
% [A, b, x, ProbInfo] = PRspherical(n, options)
% [A, b, x, ProbInfo] = PRspherical(options)
%
% This function genetates a tomography test problem based on the spherical
% Radon tranform where data consists of integrals along circles.  This type
% of problem arises, e.g., in photoacoustic imaging.
%
% The image domain is a square centered at the origin.  The centers for the
% integration circles are placed on a circle just outside the image domain.
% For each circle center we integrate along a number of concentric circles
% with equidistant radii, using the periodic trapezoidal rule.
%
% The severity of the problem increases slightly when the angular range
% decreases, and the matrix becomes singular as angles and numCircles
% become small.
%
% Input:
%  n : size of the image; a scalar such that the size is n x n.
%      Default: n = 256.
%  options : Structure containing the following optional fields:
%    phantomImage - user supplied test image of size n x n, of type numeric,
%             2-D only, or character string indicating
%             'shepplogan' : Shepp-Logan phantom
%             'smooth' : a smooth image
%             'binary' : a binary image
%             'threephases' : random image with pixel values 0, 0.5, 1
%                arranged in domains
%             'threephasessmooth' : similar to threephases, but the
%                domains have smoothly varying pixel values and there is
%                a smooth background
%             'fourphases' : similar to 'binary' but with three phases
%                separated by (thin) structures that form the fourth phase
%             'grains' : a random image with Voronoi cells
%             'ppower' : a random image with patterns of nonzero pixels
%             Default: 'shepplogan'.
%             This image is then stored in the output vector x.
%    angles     - vector of angles (in degrees) to the sources lying outside
%                 the image domain.
%                 Default is n angles between 0 and 360 degrees.
%    numCircles - number of concentric integration circles for each source.
%                 Default: numCircles = round(sqrt(2)*n).
%    sm         - logical; if true (default) then A is a sparse matrix,
%                 otherwise it is a function handle.
%
% Output:   
%  A : Sparse matrix or function handle for forward/adjoint problem.
%  b : Vector with projection data (sinogram).
%  x : Vector with image (i.e., exact image with stacked columns).
%  ProbInfo : structure containing some information about problem
%      problemType - kind of test problem (in this case: 'tomography')
%      xType       - solution type (in this case 'image2D')
%      bType       - data type (in this case 'image2D')
%      xSize       - size of image x
%      bSize       - size of sinogram b
%
% See also: PRblur, PRdiffusion, PRinvinterp2, PRnmr, PRseismic, PRtomo,
% PRnoise, PRshowb, PRshowx

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Initialization: set the default image size.
default_N = 256;

% Check for acceptable number of optional input arguments.
switch length(varargin)
    case 0
        N = default_N; options = [];
    case 1
        if isa(varargin{1}, 'double')
            N = varargin{1}; options = [];
        else
            N = default_N; options = varargin{1};
        end
    case 2
        if isa(varargin{1}, 'double')
            N = varargin{1}; options = varargin{2};
        else
            N = varargin{2}; options = varargin{1};
        end
    otherwise
        error('Too many input parameters')
end

% Set default values for options.
default_angles = linspace(0,360,N+1);
default_angles = default_angles(2:end);
defaultopt = struct('phantomImage', 'shepplogan', 'sm', true, ...
    'angles', default_angles, 'numCircles', NaN);
  
% If input is 'defaults,' return the default options in A.
if nargin == 1 && nargout <= 1 && strcmp(varargin,'defaults')
    A = defaultopt;
    return;
end

if isempty(options)
    options = defaultopt;
end

options = PRset(defaultopt, options);

phantom    = PRget(options, 'phantomImage', [], 'fast');
angles     = PRget(options, 'angles',       [], 'fast');
sm         = PRget(options, 'sm',           [], 'fast');
numCircles = PRget(options, 'numCircles',   [], 'fast');

% If a phantom image is given, then this defines N and the other parametes.
if isnumeric(phantom)
    % Make sure user input image is a matrix
    if ~ismatrix(phantom)
        error('Expected user supplied phantom image to be a 2-D array of type numeric')
    else
        N = size(phantom,1);
        if size(phantom,2)~=N, error('phantomImage image must be square'); end
        x = double(phantom(:));
        angles = linspace(0,360,N+1); angles = angles(2:end);
    end
 else
    x = phantomgallery(phantom,N);
    x = x(:);
end
if isnan(numCircles), numCircles = round(sqrt(2)*N); end

A = sphericaltomo(N,angles,numCircles,0,sm);
if sm
    b = A*x;
else
    b = A(x,'notransp');
end

% Providing information about the test problem.
ProbInfo.problemType = 'tomography';
ProbInfo.xType = 'image2D';
ProbInfo.xSize = [N N];
ProbInfo.bType = 'image2D';
ProbInfo.bSize = [numel(angles),numCircles];