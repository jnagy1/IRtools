function [A, b, x, ProbInfo] = PRblurgauss(varargin) 
% PRblurgauss Image deblurring problem with a Gaussian point spread function
%
% [A, b, x, ProbInfo] = PRblurgauss
% [A, b, x, ProbInfo] = PRblurgauss(n)
% [A, b, x, ProbInfo] = PRblurgauss(n, options)
% [A, b, x, ProbInfo] = PRblurgauss(options)
%
% Generates data for use in image deblurring problems with a Gaussian
% point spread function.
%
% Input:
%  n      -  size of the image. Can be either a scalar n (in this case 
%            the size is n x n) or a vector [nrow, ncol] (in this case the
%            size is (nrow x ncol).
%            Default: n = 256.
%  options - Structure containing the following optional fields:
%    trueImage  : test image of size n, of type numeric, 2-D only,
%                 or character string indicating
%                 'pattern1'  : geometrical image
%                 'pattern2'  : geometrical image
%                 'sppattern' : sparse (edges of ) a geometrical image
%                 'ppower'    : random image with patterns of nonzero pixels
%                 'smooth'    : very smooth image
%                 'dot2'      : two small Gaussian shaped dots, e.g., a
%                               binary star
%                 'dotk'      : n/2 small Gaussian shaped dots, e.g., stars
%                               (placement is random, reset using rng(0))
%                 'satellite' : satellite test image 
%                 'hst'       : image of the Hubble space telescope
%                 Default: 'hst'.
%                 This image is then stored in the output vector x.
%    BlurLevel  : If choosing one of the built-in PSFs, this sets the
%                 severity of the blur to one of the following:
%                 'mild'
%                 'medium'
%                 'severe'
%                 Default is 'medium'
%    BC         : Specify boundary condition:
%                 'zero'
%                 'periodic'
%                 'reflective' (or 'neumann' or 'reflexive')
%                 Default: 'reflective'
%                 Note that in this case an extended (or padded) test image
%                 is blurred using 'zero' boundary conditions, and then the 
%                 central subimage of size n is extracted from the exact and
%                 the blurred image. No inverse crime is committed, 
%                 i.e., A*x ~= b.
%    CommitCrime: To get an exact system Ax = b (i.e., commit the inverse
%                 crime), set this to:
%                 'on'
%                 Default is 'off' (do not commit the inverse crime).
%
% Output:   
%  A        - blurring matrix, which is a either a psfMatrix object in
%             the case of spatically invariant blur, or a sparse matrix
%             in the case of spatially variant blur
%  b        - blurred vector (i.e., blurred image with stacked columns)
%  x        - image vector, i.e., exact (unknown) image with stacked columns
%  ProbInfo - structure whose fields contain information about problem:
%               problemType : kind of test problem generated
%                             (in this case: 'deblurring')
%               xType       : solution type (in this case 'image2D')
%               bType       : data type (in this case 'image2D')
%               xSize       : size of image x
%               bSize       : size of image b
%               psf         : point spread function
%
% See also: PRblur, PRblurdefocus, PRblurmotion, PRblurrotation,
% PRblurshake, PRblurspeckle, PRdiffusion, PRinvinterp2, PRnmr,
% PRseismic, PRspherical, PRtomo, PRnoise, PRshowb, PRshowx, fspecial

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('trueImage', 'hst', 'BlurLevel', 'medium', ...
    'BC', 'reflective', 'CommitCrime', 'off');
  
% If input is 'defaults,' return the default options in X
if nargin == 1 && nargout <= 1 && strcmp(varargin,'defaults')
    A = defaultopt;
    return;
end

% Check for acceptable number of optional input arguments
switch length(varargin)
    case 0
        n = []; options = [];
    case 1
        if isa(varargin{1}, 'double')
            n = varargin{1}; options = [];
        else
            n = []; options = varargin{1};
        end
    case 2
        if isa(varargin{1}, 'double')
            n = varargin{1}; options = varargin{2};
        else
            n = varargin{2}; options = varargin{1};
        end
    otherwise
        error('Too many input parameters')
end

if isempty(options)
    options = defaultopt;
end

options = PRset(defaultopt, options);
options = PRset(options, 'PSF', 'gauss');
[A, b, x, ProbInfo] = PRblur(n, options);
