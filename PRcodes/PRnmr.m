function [A, b, x, ProbInfo] = PRnmr(varargin)
%PRnmr Generates data for 2D NMR relaxometry problems
%
% [A, b, x, ProbInfo] = PRnmr
% [A, b, x, ProbInfo] = PRnmr(n)
% [A, b, x, ProbInfo] = PRnmr(n, options)
%
% This function generates data for use in 2D nuclear magnetic resonance
% (NMR) relaxometry problems.
%
% NMR relaxometry is a useful tool to analyze the chemical properties of 
% porous media: given a signal sampled at different acquisition times,
% the goal of NMR is to reconstruct the joint distribution of 
% longitudinal and transverse relaxation times T1 and T2, respectively.
%
% The severity of the problem increases when Tlogleft becomes more negative;
% the other parameters have little influence.
%
% Input:
%  n       - Size of the relaxation time distribution to be recovered.
%            Must be an integer (in this case the size is n x n) or a
%            vector with two integer entries (in this case n = [n1, n2]
%            and the size is n1 x n2). Default: n = 128.
%  options - Parameters that define the problem
%              numData     : number of acquired 2D measurements; the default
%                            is m1 = 2*n1 and m2 = 2*n2
%                            [ {'double'} | m | [m1, m2] ]
%              material    : phantom for the relaxation time distribution
%                            [ {'carbonate'} | 'methane' | 'organic' | 'hydroxyl' ]
%                            This phantom is then stored in the output
%                            vector x.
%              Tloglimits  : limits for the logarithm of the relaxation times T
%                            [ {[-4, 1]} | vector with two real components ]
%                            Note: if Tloglimits = [Tlogleft, Tlogright],
%                            T1 = logspace(Tlogleft, Tlogright, n1),
%                            T2 = logspace(Tlogleft, Tlogright, n2).
%              tauloglimits: limits for the logarithm of the acquisition times tau
%                            [ {[-4, 1]} | vector with two real components ]
%                            Note: if tauloglimits = [taulogleft, taulogright],
%                            tau1 = logspace(taulogleft, taulogright, m1),
%                            tau2 = logspace(taulogleft, taulogright, m2).
%
% Output:
%  A        - Function handle for forward and adjoint problem.
%             Note: the kernel associated with 2D NMR problems is separable
%  b        - Vector with relaxation data (stacked colum-wise)
%  x        - Vector with relaxation time distribution (stacked column-wise)
%  ProbInfo - Structure containing some information about problem
%               problemType : kind of test problem generated
%                             (in this case: 'nmr')
%               xType       : solution type (in this case 'surf2D')
%               bType       : data type (in this case 'surf2D')
%               xSize       : size of image x
%               bSize       : size of image b
%               T1          : discretized T1 relaxation times
%               T2          : discretized T2 relaxation times
%               tau1        : discretized tau1 relaxation times
%               tau2        : discretized tau2 relaxation times
%               A1          : matrix with kernel acting on the T1 components
%               A2          : matrix with kernel acting on the T2 components
%
% See also: PRblur, PRdiffusion, PRinvinterp2, PRseismic, PRspherical,
% PRtomo, PRnoise, PRshowb, PRshowx

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('numData','double', 'material','carbonate',...
                    'Tloglimits',[-4, 1], 'tauloglimits',[-4, 1]);
  
% If input is 'defaults,' return the default options in X.
if nargin == 1 && nargout <= 1 && strcmp(varargin,'defaults')
    A = defaultopt;
    return;
end

% Check for acceptable number of optional input arguments
switch length(varargin)
    case 0
        dim = []; options = [];
    case 1
        if isa(varargin{1}, 'double')
            dim = varargin{1}; options = [];
        else
            dim = []; options = varargin{1};
        end
    case 2
        if isa(varargin{1}, 'double')
            dim = varargin{1}; options = varargin{2};
        else
            dim = varargin{2}; options = varargin{1};
        end
    otherwise
        error('Too many input parameters')
end

default_dim = [128, 128];

% Set a minimal size for the image.
min_dim = 16;

if isempty(dim)
    dim = default_dim;
else
    if length(dim) == 1
        dim = [dim, dim];
    elseif length(dim)>2
        error('Incorrect inputs')
    end
end

n1 = dim(1);
n2 = dim(2);

if n1 < min_dim
    n1 = min_dim;
end
if n2 < min_dim
    n2 = min_dim;
end

if isempty(options)
    options = defaultopt;
end

options       = PRset(defaultopt, options);
material      = PRget(options, 'material',      [], 'fast');
numData       = PRget(options, 'numData',       [], 'fast');
Tloglimits    = PRget(options, 'Tloglimits',    [], 'fast');
tauloglimits  = PRget(options, 'tauloglimits',  [], 'fast');

if strcmp(numData, 'double')
    m1 = 2*n1;
    m2 = 2*n2;
else
    if length(numData) == 2
        m1 = numData(1);
        m2 = numData(2);
    else
        m1 = numData;
        m2 = m1;
    end
end

if Tloglimits(1)<Tloglimits(2)
    Tlogleft  = Tloglimits(1);
    Tlogright = Tloglimits(2);
else
    Tlogleft  = Tloglimits(2);
    Tlogright = Tloglimits(1);
end

if tauloglimits(1)<tauloglimits(2)
    taulogleft  = tauloglimits(1);
    taulogright = tauloglimits(2);
else
    taulogleft  = tauloglimits(2);
    taulogright = tauloglimits(1);
end

T1log = linspace(Tlogleft, Tlogright, n1+1);
T2log = linspace(Tlogleft, Tlogright, n2+1);
% Nodes for midpoint quadrature rule (equispaced in logarithmic scale).
T1log = (T1log(1:end-1) + T1log(2:end))/2;
T2log = (T2log(1:end-1) + T2log(2:end))/2;
[T2logg, T1logg] = meshgrid(T2log, T1log);
switch material
   case {'carbonate'}
        meanT = (Tlogleft + Tlogright)/2;
        varT  = 0.1;
        x = 1/(varT^2*(2*pi))*exp(-(((T1logg - meanT)./varT).^2/2+...
            ((T2logg - meanT)./varT).^2/2));
        x(x<1e-10) = 0;
        x = x(:);
    case {'methane'}
        meanT = (Tlogleft + Tlogright)/4;
        varTa  = 0.2;
        varTb  = 0.1;
        x = 1/(varTa^2*(2*pi))*exp(-(((T1logg - meanT)./varTa).^2/2+...
            ((T2logg - 3*meanT)./varTa).^2/2))+...
            1/(varTb^2*(2*pi))*exp(-(((T1logg - 3*meanT)./varTb).^2/2+...
            ((T2logg - meanT)./varTb).^2/2));
        x(x<1e-10) = 0;
        x = x(:);
    case {'organic'}
        meanTa = (Tlogleft + Tlogright)/4;
        meanTb = (Tlogleft + Tlogright)/3;
        varTa = 0.08;
        varTb = 0.4;
        x = 1/(varTa*varTb*(2*pi))*exp(-(((T1logg - meanTa)./varTa).^2/2+...
            ((T2logg - meanTb)./varTb).^2/2));
        x(x<1e-10) = 0;
        x = x(:);       
    case {'hydroxyl'}
        meanTa = (Tlogleft + Tlogright)/4;
        meanTb = (Tlogleft + Tlogright)/6;
        varTa = 0.4;
        varTb = 0.1;
        rho = 0;
        dd = varTa^2*varTb^2 - rho^4;
        x = 1/((2*pi)*sqrt(dd))*exp(-1/2*fliplr((varTb^2/dd*(T1logg-meanTa).^2-...
            2*rho^2/dd*(T1logg-meanTa).*(T2logg-meanTb)+...
            varTa^2/dd*(T2logg-meanTb).^2)));
        x(x<1e-10) = 0;
        x = x(:);
end

tau1 = logspace(taulogleft, taulogright, m1);
tau2 = logspace(taulogleft, taulogright, m2);

[T1g, tau1g] = meshgrid(10.^T1log, tau1);
A1 = (((Tlogright-Tlogleft)/n1)*log(10)*(1 - 2*exp(-tau1g./T1g)).*T1g);
[T2g, tau2g] = meshgrid(10.^T2log, tau2);
A2 = (((Tlogright-Tlogleft)/n2)*log(10)*exp(-tau2g./T2g).*T2g);

A = @(xx,tflag) OPnmr(xx,A1,A2,tflag);

b = A(x, 'notransp');

ProbInfo.problemType = 'nmr';
ProbInfo.xType = 'surf2D';
ProbInfo.bType = 'surf2D';
ProbInfo.xSize = [n1,n2];
ProbInfo.bSize = [m1,m2];
ProbInfo.T1 = 10.^T1log;
ProbInfo.T2 = 10.^T2log;
ProbInfo.tau1 = tau1;
ProbInfo.tau2 = tau2;
ProbInfo.A1 = A1;
ProbInfo.A2 = A2;