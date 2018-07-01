function [A, b, x, ProbInfo] = PRseismic(varargin) 
% PRseismic Generates data for seismic travel-time tomography problems
%
% [A, b, x, ProbInfo] = PRseismic
% [A, b, x, ProbInfo] = PRseismic(n)
% [A, b, x, ProbInfo] = PRseismic(n, options)
% [A, b, x, ProbInfo] = PRseismic(options)
%
% This function creates a 2D seismic travel-time tomography test problem
% with an N-times-N pixel domain, using s sources located on the right
% boundary and p receivers (seismographs) scattered along the left and top
% boundaries.  The rays are transmitted from each source to each receiver.
%
% The severity of the problem is almost unaffected by the problem parameters.
%
% Input:
%  n : size of the image; a scalar such that the size is n x n.
%      Default: N = 256.
%  options : Structure containing the following optional fields:
%    phantomImage - user supplied test image of size n x n, of type numeric,
%             2-D only, or character string indicating
%             'tectonic' : two tectonic plates create a subduction zone
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
%             Default: 'tectonic'.
%             This image is then stored in the output vector x.
%    wavemodel - string that defines the type of problem:
%                  'ray'     : waves are represented by straight lines (default),
%                  'fresnel' : waves travel within the first Fresnel zone.
%    s     - number of sources in the right side of the domain. Default s = n.
%    p     - number of receivers (seismographs) equally spaced on the surface
%            and on the left side of the domain. Default p = 2*n.
%    omega - Fresnel model only: dominant frequency of the propagating wave.
%            Default = 10.
%    sm    - logical; if true (default) then A is a sparse matrix, otherwise
%            it is a function handle.
%
% Output:   
%  A : Sparse matrix or function handle for forward/adjoint problem.
%  b : Vector with projection data (sinogram).
%  x : Vector with image (i.e., exact image with stacked columns).
%  ProbInfo : structure containing some information about problem
%      problemType - kind of test problem (in this case: 'tomography')
%      xType - solution type (in this case 'image2D')
%      bType - data type (in this case 'image2D')
%      xSize - size of image x
%      bSize - size of sinogram b
%
% See also: PRblur, PRdiffusion, PRinvinterp2, PRnmr, PRspherical, PRtomo,
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
defaultopt = struct('phantomImage', 'tectonic', 'sm', true, ...
    'wavemodel', 'ray', 's', N, 'p', NaN, 'omega', 10);
  
% If input is 'defaults,' return the default options in A.
if nargin == 1 && nargout <= 1 && strcmp(varargin,'defaults')
    A = defaultopt;
    return;
end

if isempty(options)
    options = defaultopt;
end

options = PRset(defaultopt, options);

phantom   = PRget(options, 'phantomImage', [], 'fast');
sm        = PRget(options, 'sm',           [], 'fast');
wavemodel = PRget(options, 'wavemodel',    [], 'fast');
s         = PRget(options, 's',            [], 'fast');
p         = PRget(options, 'p',            [], 'fast');
omega     = PRget(options, 'omega',        [], 'fast');

% If a phantom image is given, then this defines N and the other parametes.
if isnumeric(phantom)
    % Make sure user suppled phantom image is a matrix
    if ~ismatrix(phantom)
        error('Expected user supplied phantom image to be a 2-D array of type numeric')
    else
        N = size(phantom,1);
        if size(phantom,2)~=N, error('phantomImage image must be square'); end
        x = double(phantom(:));
        s = N;
    end
 else
    x = phantomgallery(phantom,N);
    x = x(:);
end
if isnan(p), p = 2*N; end

switch wavemodel
    case 'ray'
        A = seismictomo(N,s,p,0,sm);
    case 'fresnel'
        A = seismicwavetomo(N,s,p,omega,0,sm);
   otherwise
        error('Type of seismic problem not provided')
end

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
ProbInfo.bSize = [p,s];