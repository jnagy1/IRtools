function [A, b, x, ProbInfo] = PRdiffusion(varargin) 
% PRdiffusion Generates data for use in an inverse diffusion problems
%
% [A, b, x, ProbInfo] = PRdiffusion
% [A, b, x, ProbInfo] = PRdiffusion(n)
% [A, b, x, ProbInfo] = PRdiffusion(n, options)
%
% This function generates a 2D inverse diffusion problem on a uniform
% finite-element mesh with 2*(n-1)^2 triangular elements; think of the
% domain as an n-by-n pixel grid with two triangular elements in each pixel.
%
% The data represents a function u(t) on the grid at time t = Tfinal.
% The goal of the inverse problem is to compute the function u(t) at
% time t = 0 which, when diffused, produces the data.
%
% The forward problem is a very simple diffusion problem
%    du/dt = u_xx + u_yy,   0 < t < Tfinal
% with homogenous Neumann boundary conditions.  It is solved by means of
% the Crank-Nicolson-Galerkin finite-element method.  The inverse problem
% is to compute the initial condition, the function u(0) at time t = 0,
% from the solution u(t) at time t = Tfinal.
%
% The functions u(0) and u(Tfinal) are represented by the vecors x and b,
% resp., which hold the values at the corners of the finite elements.
% Both vectors have a total of n^2 elements.
%
% The severity of the problem increases with Tfinal and does not depend on
% Tsteps.
%
% Input:
%  n       - The grid has n-1 pixels in each direction, and n must be an
%            integer.  Default: n = 128.  There are 2*(n-1)^2 finite
%            elements, and the number of unknowns is n^2.
%  options - Parameters that define the problem
%              Tfinal : diffusion time; default Tfinal = 0.01.
%              Tsteps : number of time steps; default Tsteps = 100.
%
% Output:
%  A -  Function handle for the forward and adjoint problem.
%  b -  Vector that represents the solution u(t) at t = Tfinal.
%  x -  Vector that represents the initial condition u(t) at t = 0.
%  ProbInfo -  Structure containing some information about problem
%                problemType : type of test problem generated
%                              (in this case: 'diffusion')
%                xType       : solution type (in this case 'fem')
%                bType       : data type (in this case 'fem')
%                elmtab      : table of FEM elements
%                elmX        : X-coordinates of element grid
%                elmY        : Y-coordinates of element grid
%
% See also: PRblur, PRinvinterp2, PRnmr, PRseismic, PRspherical, PRtomo,
% PRnoise, PRshowb, PRshowx

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% Based on original code by Allan P. Engsig-Karup, DTU Compute.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('Tfinal', 0.01, 'Tsteps', 100);
  
% If input is 'defaults,' return the default options in A.
if nargin == 1 && nargout <= 1 && strcmp(varargin,'defaults')
    A = defaultopt;
    return;
end

% Check for acceptable number of optional input arguments.
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
Tfinal = PRget(options, 'Tfinal', [], 'fast');
Tsteps = PRget(options, 'Tsteps', [], 'fast');

if isempty(n), n = 128; end

% Create rectangular mesh consisting of uniformly distributed triangles.
[elmX,elmY] = meshgrid(linspace(0,1,n));
elmX = elmX(:);
elmY = elmY(:);

% Element table for uniformly distributed triangles in a rectangular domain.
elmtab = zeros(2*(n-1)^2,3);
count = 0;
for j = 1:n-1
    for i = 1:n-1
        elmtab(count+(1:2),:) = [i   i+1+n i+n;
                                 i+1 i+1+n i     ] + (j-1)*n;
        count = count+2;
    end
end

% Initial condition.
x = 0.7*exp( -((elmX-0.4)/0.12).^2 - ((elmY-0.5)/0.15).^2) + ...
        exp( -((elmX-0.7)/0.1 ).^2 - ((elmY-0.4)/0.08).^2);

% Create the matrices.
dt = Tfinal/Tsteps;
[AA,CC] = assembly(elmX,elmY,elmtab);
S = CC - 0.5*dt*AA;
R = CC + 0.5*dt*AA;

% A is a function handle to the function that carries out the forward
% operation, i.e., interpolates from the regular grid to the randon points,
% or the corresponding adjoint operation (representing A').
A = @(xx,tflag) OPdiffusion(xx,R,S,Tsteps,tflag);

% The right-hand side (the data) ...
b = A(x,'notransp');

ProbInfo.problemType = 'diffusion';
ProbInfo.bType = 'fem';
ProbInfo.xType = 'fem';
ProbInfo.elmtab = elmtab;
ProbInfo.elmX = elmX;
ProbInfo.elmY = elmY;
ProbInfo.R = R;                  %PCH
ProbInfo.S = S;

% SUBFUNCTIONS -----------------------------------------------------------

function [A,C] = assembly(x,y,elmtab)
% Generate linear system coefficient matrix and right hand side vector.
% Symmetry is not exploited.

P = size(elmtab,1);
iL = zeros(3*3*P,1);
jL = iL;
AL = iL;
CL = iL;
count = 0;
for p = 1:P
    [delta,abc] = basfun(p,x,y,elmtab);
    for r = 1:3
        i = elmtab(p,r);
        for s = 1:3   % Does not exploit symmetry.
            j = elmtab(p,s);
            Ars = 1/(4*abs(delta))*(abc(r,2)*abc(s,2) + abc(r,3)*abc(s,3));
            if r == s
                Crs = abs(delta)/6;
            else
                Crs = abs(delta)/12;
            end
            count = count + 1;
            iL(count) = i;
            jL(count) = j;
            AL(count) = Ars;
            CL(count) = Crs;
        end
    end        
end
iL = iL(1:count);
jL = jL(1:count);
AL = AL(1:count);
A = sparse(iL,jL,AL);
C = sparse(iL,jL,CL);

function [delta,abc] = basfun(p,x,y,elmtab)
% For a given element compute the geometric properties of the element.

xc = x(elmtab(p,:));
yc = y(elmtab(p,:));
delta = 0.5*( xc(2)*yc(3) - yc(2)*xc(3) - (xc(1)*yc(3) - ...
              yc(1)*xc(3)) + xc(1)*yc(2) - yc(1)*xc(2)      );

% Vectorized output.
j   = [2 3 1]';
k   = [3 1 2]';
at  = xc(j).*yc(k) - xc(k).*yc(j);
bt  = yc(j) - yc(k);
ct  = xc(k) - xc(j);
abc = [at bt ct];