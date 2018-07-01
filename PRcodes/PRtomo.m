function [A, b, x, ProbInfo] = PRtomo(varargin) 
% IRtomo Generates data for X-ray tomographic reconstruction problems
%
% [A, b, x, ProbInfo] = PRtomo
% [A, b, x, ProbInfo] = PRtomo(n)
% [A, b, x, ProbInfo] = PRtomo(n, options)
% [A, b, x, ProbInfo] = PRtomo(options)
%
% This function uses the "line model" to create a 2D X-ray tomography test
% problem with an N-times-N pixel domain, using p rays for each angle in
% the vector theta.
%
% The severity of the problem is almost unaffected by the problem parameters.
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
%    CTtype - string that defines the type of CT problem:
%             'parallel'  : parallel beam geometry (default),
%             'fancurved' : fan beam with curved detector.
%    sm     - logical; if true (default )then A is a sparse matrix,
%             otherwise it is a function handle.
%    angles - vector of projection angles, in degrees.
%             Default: 0:1:179 for parallel beam, 0:2:358 for fan beam.
%    p      - number of rays for each source angle.
%             Default: p = round(sqrt(2)*n).
%    d      - Parallel beam only: scalar denoting the distance from the
%             first ray to the last; default: d = p-1.
%    R      - Fan beam only: the distance from the source to the center
%             of the phantom is R*N.  Default: R = 2.
%    span   - Fan beam only: scalar that determines the angular span of
%             the rays, in degrees. The default value is defined such
%             that from the source at (0,R*N) the first and last rays
%             hit the image corners (-n/2,n/2) and (n/2,n/2).
%
% Output:   
%  A : Sparse matrix or function handle for forward/adjoint problem.
%  b : Vector with projection data (the sinogram).
%  x : Vector with image (i.e., exact image with stacked columns).
%  ProbInfo : structure containing some information about problem
%      problemType - kind of test problem (in this case: 'tomography')
%      xType       - solution type (in this case 'image2D')
%      bType       - data type (in this case 'image2D')
%      xSize       - size of image x
%      bSize       - size of sinogram b
%
% See also: PRblur, PRdiffusion, PRinvinterp2, PRnmr, PRseismic,
% PRspherical, PRtomo, PRnoise, PRshowb, PRshowx, radon

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
defaultopt = struct('phantomImage', 'shepplogan', 'CTtype', 'parallel', ...
    'sm', true, 'angles', 0:1:179, 'p', NaN, 'R', 2, 'd', NaN, 'span', NaN);
  
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
CTtype     = PRget(options, 'CTtype',       [], 'fast');
angles     = PRget(options, 'angles',       [], 'fast');
sm         = PRget(options, 'sm',           [], 'fast');
p          = PRget(options, 'p',            [], 'fast');
R          = PRget(options, 'R',            [], 'fast');
d          = PRget(options, 'd',            [], 'fast');
span       = PRget(options, 'span',         [], 'fast');

% If a phantom image is given, then this defines N and the other parametes.
if isnumeric(phantom)
    % Make sure user input image is a matrix
    if ~ismatrix(phantom)
        error('Expected user supplied phantom image to be a 2-D array of type numeric')
    else
        N = size(phantom,1);
        if size(phantom,2)~=N, error('phantomImage image must be square'); end
        x = double(phantom(:));
    end 
else
    x = phantomgallery(phantom,N);
    x = x(:);
end
if isnan(p), p = round(sqrt(2)*N); end

switch CTtype
    case 'parallel'
        if isnan(d), d = p-1; end
        A = paralleltomo(N,angles,p,d,0,sm);
    case 'fancurved'
        if isnan(span), span = 2*atand(1/(2*R-1)); end
        A = fancurvedtomo(N,angles,p,R,span,0,sm);
   otherwise
        error('Type of CT problem not provided')
end
if sm
    b = A*x;
else
    b = A(x,'notransp');
end

% Providing information about the test problem.
ProbInfo.problemType = 'tomography';
ProbInfo.xType = 'image2D';
ProbInfo.bType = 'image2D';
ProbInfo.xSize = [N,N];
ProbInfo.bSize = [p,length(angles)];