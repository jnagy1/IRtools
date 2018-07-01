function [A, b, x, ProbInfo] = PRinvinterp2(varargin) 
% PRinvinterp2 Generates data for use in 2D inverse interpolation problems
%
% [A, b, x, ProbInfo] = PRinvinterp2
% [A, b, x, ProbInfo] = PRinvinterp2(n)
% [A, b, x, ProbInfo] = PRinvinterp2(n, options)
%
% This function generates a 2D inverse interpolation problem on an n-by-n
% regular grid.  Given n^2 data on random points the goal is to compute
% function values of the grid which, when interpolated, produces the data.
%
% The decay of the singular values is unaffected by the interpolation method,
% but the matrix is rank deficient for 'nearest' and 'linear' interpolation.
%
% Input:
%  n       - The size of the image is n x n, and n must be an integer.
%            Default: n = 128.
%  options - Parameters that define the problem
%              InterpMethod - String that defines the interpolation:
%              'nearest', 'linear' (default), 'cubic' or 'spline'.
%
% Output:
%  A        - Function handle for forward and adjoint problem.
%  b        - Vector with function values at random points.
%  x        - Vector with function values on regular mesh.
%  ProbInfo - Structure containing some information about problem
%               problemType : kind of test problem generated
%                             (in this case: 'invinterp2')
%               xType       : solution type (in this case 'surf2D')
%               bType       : data type (in this case 'surf2D')
%               xSize       : size of image x
%               bSize       : size of image b
%               Xr          : x-coordinates of data b
%               Yr          : y-coordinates of data b
%
% Note that A(x,'notransp') does not equal b because of rounding errors.
%
% See also: PRblur, PRdiffusion, PRnmr, PRseismic, PRspherical, PRtomo,
% PRnoise, PRshowb, PRshowx

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('InterpMethod', 'linear');
  
% If input is 'defaults,' return the default options in A.
if nargin == 1 && nargout <= 1 && strcmp(varargin,'defaults')
    A = defaultopt;
    return;
end

% Check for acceptable number of optional input arguments
switch length(varargin)
    case 0
        N = []; options = [];
    case 1
        if isa(varargin{1}, 'double')
            N = varargin{1}; options = [];
        else
            N = []; options = varargin{1};
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

if isempty(options)
    options = defaultopt;
end

options = PRset(defaultopt, options);
method  = PRget(options, 'InterpMethod', [], 'fast');

if isempty(N), N = 128; end

% Define the function to be interpolated.
f = @(x,y) sin(pi*x) .* sin(0.5*pi*y);

% The exact solution consists of the function values af regular grid points.
X = linspace(0,1,N);
Y = linspace(0,1,N);
[XX,YY] = meshgrid(X,Y);
x = f(XX,YY);
x = x(:);

% The right-hand side (the data) consists of the function values at
% random points.
Xr = rand(N^2,1);
Yr = rand(N^2,1);
b = f(Xr,Yr);
b = b(:);

% A is a function handle to the function that carries out the forward
% operation, i.e., interpolates from the regular grid to the randon points,
% or the corresponding adjoint operation (representing A').
A = @(xx,tflag) OPinvinterp2(xx,X,Y,Xr,Yr,method,tflag);

ProbInfo.problemType = 'invinterp2';
ProbInfo.bType = 'surf2D';
ProbInfo.xType = 'surf2D';
ProbInfo.xSize = [N,N];
ProbInfo.bSize = [N,N];
ProbInfo.Xr = Xr;
ProbInfo.Yr = Yr;