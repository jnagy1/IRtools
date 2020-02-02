function [bn, NoiseInfo] = PRnoise(b, varargin)
% PRnoise Add noise of a specified distribution and level
%
% [bn, NoiseInfo] = PRnoise(b);
% [bn, NoiseInfo] = PRnoise(b, kind);
% [bn, NoiseInfo] = PRnoise(b, NoiseLevel);
% [bn, NoiseInfo] = PRnoise(b, kind, NoiseLevel);
% [bn, NoiseInfo] = PRnoise(b, NoiseLevel, kind);
%
% This function adds noise of a specified kind and level to data, such that
% bn = b + noise and ||noise(:)||_2 / ||b(:)||_2 = NoiseLevel .
%
% Input:
%   b          - data vector or array
%   kind       - kind of noise
%                [ {'gauss'} | 'laplace' | 'multiplicative' ]
%   NoiseLevel - relative level of noise, defined as ||noise(:)||/||b(:)||
%                [ {0.01} | nonnegative scalar ]
%
% Output:
%   bn        - the noisy data bn = b + noise
%   NoiseInfo - structure whose fileds contain information about the noise
%                 kind       : kind of test problem generated
%                 NoiseLevel : relative level of noise
%                 noise      : the noise vector or array
%
% Gaussian noise is created by means of MATLAB's randn function.
%
% Laplacian noise is generated as follows:
%    r = rand(numel(b),1);
%    r = sign(0.5-r).*(1/sqrt(2)).*log(2*min(r,1-r));
%    noise = ((NoiseLevel*norm(b(:)))/norm(r))*r;
%
% For multiplicative noise, each element of bn equals the corresponding
% element of b times a random variable following a Gamma distribution
% with mean 1, scaled to approximately obtain the desired noise level.
% The Statistics Toolbox is required for this type of noise.
%
% Other types of noise may be added with the function imnoise from the
% Image Processing Toolbox.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Initialization: set the default inputs
default_kind = 'gauss';
default_level = 0.01;

if nargin == 0
    error('Not enough input arguments')
elseif nargin > 3
    error('Too many input arguments')
end

switch length(varargin)
    case 0
        kind = []; NoiseLevel = [];
    case 1
        if isa(varargin{1}, 'double')
            kind = []; NoiseLevel = varargin{1};
        else
            kind = varargin{1}; NoiseLevel = [];
        end
    case 2
        if isa(varargin{1}, 'double')
            kind = varargin{2}; NoiseLevel = varargin{1};
        else
            kind = varargin{1}; NoiseLevel = varargin{2};
        end
end

if isempty(kind)
    kind = default_kind;
end
if isempty(NoiseLevel)
    NoiseLevel = default_level;
end

if NoiseLevel<0, error('The level of noise must be nonnegative'), end

[m,n] = size(b); 
M = m*n;

if strcmpi(kind, 'gauss')
    r = randn(M,1);
    noise = ((NoiseLevel*norm(b(:)))/norm(r))*r;
    noise = reshape(noise,m,n);
    bn = b + noise;
elseif strcmpi(kind, 'laplace')
    r = rand(M,1);
    r = sign(0.5-r).*(1/sqrt(2)).*log(2*min(r,1-r));
    noise = ((NoiseLevel*norm(b(:)))/norm(r))*r;
    noise = reshape(noise,m,n);
    bn = b + noise;
elseif strcmpi(kind, 'multiplicative')
    if ~license('test','Statistics_toolbox')
        error('It appears that the Statistics Toolbox is not available')
    end
    if NoiseLevel==0
        bn = b;
    else
        dfun = @(k) norm(b(:).*mygamrnd(k,1/k,M)-b(:))/norm(b(:))-NoiseLevel;
        options.TolX = 1e-4;
        kappa = fzero(dfun,1/NoiseLevel^2,options);
        bn = b.*reshape(mygamrnd(kappa,1/kappa,M),m,n);
    end
    noise = bn - b;
else
    error('Invalid noise type')
end

NoiseInfo.kind  = kind;
NoiseInfo.NoiseLevel = NoiseLevel;
NoiseInfo.noise = noise;

% Subfunction ===========================================================

function Y = mygamrnd(a,b,N)
%mygamrnd  Gamma random numbers (via tranformation-rejection)
%
%  Y = mygamrnd(a,b,N)
%
% Produces a column vector of length N with random numbers for a Gamma
% distribution with shape parameter A and scale parameter B.
% This is an implementation of the transformation-rejection method:
% G. Marsaglia and W. W. Tsang, A simple method for generating Gamma
% variables, ACM Trans. Math. Soft., 26 (2000), pp. 363–372.

% The case a < 1 is handled as a special case.
if a < 1
    Y = mygamrnd(a+1,b,N).*(rand(N,1).^(1/a));
    return
end

% Need to iterate since each pass is not guaranteed to produce N values
% satisfying the acceptance criterion.
Y = zeros(0,1);
while length(Y) < N
    d = a - 1/3;
    X = randn(N,1);
    U = rand(N,1);
    V = (1+X./sqrt(9*d)).^3;
    accept = (V > 0) & log(U) < (0.5*X.^2 + d - d*V + d*log(V));
    YY = d*(V(accept)).*b;
    Y = [Y;YY];
end

% Return a vector of the desired length.
Y = Y(1:N);