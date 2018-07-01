function [X,info] = IRart(A,b,varargin)
%IRart Algebraic Reconstruction Technique
%
% options  = IRart('defaults')
% [X,info] = IRart(A,b,K)
% [X,info] = IRart(A,b,K,options)
%
% This function calls the function 'art' in AIR Tools II which implements
% the Algebraic Reconstruction Technique (ART) method, also known as
% Kaczmarz's method.
%
% With 'defaults' as input returns the default options.  Otherwise outputs
% the iterates specified in K, using max(K) as MaxIter, and using all other
% default options.  With options as input: uses the user-specified options
% and all the other default options.
%
% Inputs:
%  A : either (a) a full or sparse matrix
%             (b) a matrix object that performs the matrix*vector operation
%             (c) user-defined function handle
%  b : right-hand side vector
%  K : (optional) integer vector that specifies which iterates are returned
%      in X; the maximum number of iterations is assumed to be max(K)
%      [ positive integer | vector of positive components ]
%  options : structure with the following fields (optional)
%      x0            - initial guess for iterations; default = zero vector
%                      [ array | {'none'} ]
%      MaxIter       - maximum allowed number of iterations
%                      [ positive integer | {100} ]
%      x_true        - true solution; allows us to returns error norms
%                      with respect to x_true at each iteration
%                      [ array | {'none'} ]
%      stopCrit      - stopping criterion for the iterations
%                      [ {'none'} | 'discrep' ]
%                      Note: 'discrepancy' requires NoiseLevel and eta
%      NoiseLevel    - norm of noise in rhs divided by norm of rhs 
%                      [ {'none'} | nonnegative scalar ]
%      eta           - safety factor for the discrepancy principle
%                      [ {1.01} | scalar greater than (and close to) 1 ]
%      relaxParam    - constant relaxation parameter bounded above by 2
%                      [ positive scalar | {'none'} ]
%                      If 'none' then a good value is chosed by the function
%      nonnegativity - apply nonnegativity constraints
%                      [ 'on' | {'off'} ]
%      Ubound        - apply box constraints in the interval [0,Ubound]
%                      [ positive scalar | {'off'} ]
%      IterBar       - shows the progress of the iterations
%                      [ {'on'} | 'off' ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X    : computed solutions, stored column-wise at the iterations listed in K
%   info : structure with the following fields:
%      its              - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag         - a flag that describes the stopping condition:
%                         1 : reached maximum number of iterations
%                         2 : discrepancy principle satisfied
%      relaxParam       - the used relaxation parameter
%      Enrm             - relative error norms (requires x_true) for the
%                         stored iterations
%
% Note that this function provides a simplified call to the function "art"
% in AIR Tools II which has more features than we allow here.  To use the
% full power of the ART methods consider using AIR Tools II directly.
%
% See also: IRsirt, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0','none', 'IterBar','on', 'stopCrit','none', ...
    'NoiseLevel','none', 'eta',1.01, 'nonnegativity','off', ...
    'Ubound','none', 'relaxParam','none', 'MaxIter',100, 'x_true','none');
  
% If input is 'defaults,' return the default options in X.
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    X = defaultopt;
    return
end

% Check for acceptable number of optional input arguments.
switch length(varargin)
    case 0
        K = []; options = [];
    case 1
        if isa(varargin{1}, 'double')
            K = varargin{1}; options = [];
        else
            K = []; options = varargin{1};
        end
    case 2
        if isa(varargin{1}, 'double')
            K = varargin{1}; options = varargin{2};
        else
            K = varargin{2}; options = varargin{1};
        end
    otherwise
        error('Too many input parameters')
end

if isempty(options)
    options = defaultopt;
end

options = IRset(defaultopt, options);

if min(K) < 1, error('Number of iterations must be positive'), end

x0 = IRget(options, 'x0', [], 'fast');
if strcmpi(x0, 'none'), x0 = []; end

MaxIter       = IRget(options, 'MaxIter',       [], 'fast');
IterBar       = IRget(options, 'IterBar',       [], 'fast');
relaxParam    = IRget(options, 'relaxParam',    [], 'fast');
stopCrit      = IRget(options, 'stopCrit',      [], 'fast');
nonnegativity = IRget(options, 'nonnegativity', [], 'fast');
Ubound        = IRget(options, 'Ubound',        [], 'fast');
x_true        = IRget(options, 'x_true',        [], 'fast');

if isempty(K)
    K = MaxIter;
end
% Sorting the iteration numbers (in case they are shuffled in input).
K = sort(K,'ascend'); K = unique(K);
if ~((isreal(K) && (all(K > 0)) && all(K == floor(K))))
    error('K must be a vector of positive real integers')
end

% Set options for ART.  Always used a fixed lambda, either chosen by the
% user or set by the ART function.  Always use no stopping rule or stop
% by the size of the residual (discrepancy principle).

if strcmpi(relaxParam,'none')
    % The field lambda is not present.
else
   artoptions.relaxpar = relaxParam;
end

if strcmpi(stopCrit,'discrep')
    artoptions.stoprule.type = 'DP';
    if isempty(options.NoiseLevel) || strcmp(options.NoiseLevel,'none')
        error('When using discrepancy principle: options.NoiseLevel must be specified')
    end
    artoptions.stoprule.taudelta = options.eta*options.NoiseLevel*norm(b);
else
    artoptions.stoprule.type = 'none';
end

if strcmpi(nonnegativity,'on')
    artoptions.lbound = 0;
end

if (isreal(Ubound) && isscalar(Ubound) && Ubound > 0)
    artoptions.lbound = 0;
    artoptions.ubound = Ubound;
end

if strcmp(IterBar,'on')
    artoptions.waitbar = true;
end
    
% Call the ART method.
[X,artinfo] = kaczmarz(A,b,K,x0,artoptions);

switch artinfo.stoprule
    case 0
        info.StopFlag = 1;
    case 2
        info.StopFlag = 2;
        if artinfo.finaliter==0
            warning(['No iterations were performed, beause the starting',...
                     ' vector x0 satisfies the discrepancy principle'])
        end
end
info.its = artinfo.finaliter;
info.saved_iterations = artinfo.itersaved;
info.relaxParam = artinfo.relaxpar;

% Compute relative error norms, if requested.
if ~strcmp(x_true,'none')
    Enrm = zeros(length(info.saved_iterations),1);
    for k=1:length(info.saved_iterations)
        Enrm(k) = norm(x_true-X(:,k));
    end
    info.Enrm = Enrm/norm(x_true);
end