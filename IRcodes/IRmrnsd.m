function [X, info] = IRmrnsd(A, b, varargin)
%IRmrnsd Modified Residual Norm Steepest Descent method
%
% options  = IRmrnsd('defaults')
% [X,info] = IRmrnsd(A,b)
% [X,info] = IRmrnsd(A,b,K)
% [X,info] = IRmrnsd(A,b,options)
% [X,info] = IRmrnsd(A,b,K,options)
%
% This function implements the Modified Residual Norm Steepest Descent
% method for computing a nonnegatiely constrained least squares solution.
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
%      x0         - Initial guess for the iterations
%                   [ array | {'none'} ]
%      MaxIter    - maximum allowed number of cgls iterations
%                   [ {100 } | positive integer ]
%                   NOTE: K overrules MaxIter if both are assigned
%      x_true     - true solution; allows us to returns error norms with
%                   respect to x_true at each iteration
%                   [ array | {'none'} ]
%      NoiseLevel - norm of noise in rhs divided by norm of rhs 
%                   [ {'none'} | nonnegative scalar]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      NE_Rtol    - relative tolerance on the normal equation residual norm
%                   [ {1e-12} | positive integer ]
%      NoStop     - specifies whether the iterations should proceed
%                   after a stopping criterion has been satisfied
%                   [ 'on' | {'off'} ]
%      IterBar    - shows the progress of the iterations
%                   [ {'on'} | 'off' ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its      - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - string that describes the inner stopping condition:
%                   * Residual tolerance satisfied (discrepancy principle)
%                   * Normal equations residual tolerance satisfied
%                   * Reached maximum number of iterations
%      Rnrm     - relative residual norms at each iteration
%      NE_Rnrm  - normal eqs relative residual norms at each iteration
%      Xnrm     - solution norms at each iteration
%      Enrm     - relative error norms (requires x_true) at each iteration
%      StopReg  - struct containing information about the solution that
%                 satisfies the stopping criterion.  Fields:
%                   It   : iteration where the stopping criterion is satisfied
%                   X    : solution satisfying the stopping criterion
%                   Enrm : the corresponding relative error (requires x_true)
%      BestReg  - struct containing information about the solution that
%                 minimizes Enrm (requires x_true), with the fields:
%                   It   : iteration where the minimum is attained
%                   X    : best solution
%                   Enrm : best relative error
%
% See also: IRconstr_ls, IRfista, IRnnfcgls, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0','none', 'MaxIter',100, 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, 'IterBar','on', ...
    'NoStop', 'off');
  
% If input is 'defaults,' return the default options in X.
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    X = defaultopt;
    return;
end

defaultopt.verbosity = 'on';

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
        if isfield(options, 'MaxIter') && ~isempty(options.MaxIter) && (~isempty(K) && options.MaxIter ~= max(K))
            warning('The value of MaxIter is discarded; the maximum value in K is taken as MaxIter')
        end
    otherwise
        error('Too many input parameters')
end

if isempty(options)
    options = defaultopt;
end

options = IRset(defaultopt, options);

MaxIter    = IRget(options, 'MaxIter',    [], 'fast');
x_true     = IRget(options, 'x_true',     [], 'fast');
NoiseLevel = IRget(options, 'NoiseLevel', [], 'fast');
eta        = IRget(options, 'eta',        [], 'fast');
NE_Rtol    = IRget(options, 'NE_Rtol',    [], 'fast');
IterBar    = IRget(options, 'IterBar',    [], 'fast');
NoStop     = IRget(options, 'NoStop',     [], 'fast');
verbose    = IRget(options, 'verbosity',  [], 'fast');

verbose = strcmp(verbose, 'on');

% Setting K.
if isempty(K)
    K = MaxIter;
end
% Sorting the iterations (in case they are shuffled in input).
K = K(:); K = sort(K,'ascend'); K = unique(K);
if ~((isreal(K) && (all(K > 0)) && all(K == floor(K))))
    error('K must be a vector of positive real integers')
end
if K(end) ~= MaxIter
    MaxIter = K(end);  
end

StopIt = MaxIter;

if isempty(NoiseLevel) || strcmp(NoiseLevel,'none')
    Rtol = 0;
else
    Rtol = eta*NoiseLevel;
end

% We need to find the number of columns in matrix A, but if A is not given 
% as a matrix, and no initial guess is given, then we can find it by 
% computing A'*b.  Since we need this anyway, it doesn't cost any 
% additional work.
trAb = Atransp_times_vec(A, b);
n = length(trAb);

nrmb = norm(b(:));
nrmAtb = norm(trAb(:));

% See if an initial guess is given.  If not, then use 0 as initial guess.  
x = IRget(options, 'x0', [], 'fast');

if strcmp(x,'none')
    % the default initial guess for the iterations is defined as the
    % minimizer of || b - A*(alpha*ones(n,1)) ||_2
    coeffx0 = A_times_vec(A, ones(n,1)); alpha = (coeffx0'*b)/norm(coeffx0)^2;
    if alpha <= 0
        alpha = sqrt(eps);
    end
    x = alpha*ones(n,1);
end

if norm(x) == 0
    error(['IRmrnsd cannot handle a zero initial guess. ',...
    'Please consider assigning a different initial guess, or avoid specifying the initial guess in order to have a default value'])
end

% If initial guess has negative values, compensate.
minx = min(x(:));
if minx < 0
     x = x - min(0,minx) + sqrt(eps);
end

noIterBar = strcmp(IterBar,{'off'});

% Declare matrices.
X = zeros(n,length(K));
Xnrm    = zeros(MaxIter,1);
Rnrm    = zeros(MaxIter,1);
NE_Rnrm = zeros(MaxIter,1);
if strcmp(x_true,'none')
    errornorms = false;
else
    errornorms = true;
    Enrm = zeros(MaxIter,1);
    nrmtrue = norm(x_true(:));
    BestReg.It = [];
    BestReg.X = [];
    BestReg.Enrm = [];
    BestEnrm = 1e10;
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestReg.NE_Rnrm = [];
end

NoStop = strcmp(NoStop,'on');
saved_iterations = zeros(1, length(K));

% Initializing some variables.
r = b - A_times_vec(A,x);
g = -Atransp_times_vec(A, r);
xg = x .* g;
gamma = g(:)' * xg(:);

if ~noIterBar
  h_wait = waitbar(0, 'Running iterations, please wait ...');
end

j = 0;
for k = 1:MaxIter
    if ~noIterBar
        waitbar(k/MaxIter, h_wait)
    end
  
    s = - x .* g;

    u = A_times_vec(A, s);
  
    theta = gamma / (u(:)'*u(:));
    neg_ind = s < 0;
  
    alpha = min( theta, min( -x(neg_ind) ./ s(neg_ind) ) );
    if isempty(alpha)
        alpha = theta;
    end
  
    x = x + alpha*s;
    
    r = r - alpha*u;
    AlreadySaved = 0; 
    if any(K == k)
        j = j+1;
        X(:,j) = x;
        saved_iterations(j) = k;
        AlreadySaved = 1; 
    end
    
    z = Atransp_times_vec(A, u);

    g = g + alpha*z;
    xg = x .* g;
    gamma = g(:)' * xg(:);

    % Compute norms.
    Xnrm(k)    = norm(x(:));
    Rnrm(k)    = norm(r(:))/nrmb;
    NE_Rnrm(k) = sqrt(gamma)/nrmAtb;
    if errornorms
        Enrm(k) = norm(x_true-x)/nrmtrue;
        if Enrm(k)<BestEnrm
            BestReg.It = k;
            BestReg.X = x;
            BestEnrm = Enrm(k);
            BestReg.Enrm = BestEnrm;
            BestReg.Xnrm = Xnrm(k);
            BestReg.Rnrm = Rnrm(k);
            BestReg.NE_Rnrm = NE_Rnrm(k);
        end
    end  
    if (Rnrm(k) <= Rtol)  && (StopIt == MaxIter)
        if verbose
            disp('Residual tolerance satisfied')
        end
        StopFlag = 'Residual tolerance satisfied';
        StopIt = k;
        StopReg.It = k;
        StopReg.X = x;
        if errornorms, StopReg.Enrm = Enrm(k); end
        if ~ NoStop
            if ~AlreadySaved
                j = j+1;
                X(:,j) = x;
                saved_iterations(j) = k;
                AlreadySaved = 1;
            end
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if errornorms, Enrm = Enrm(1:k); end
            X = X(:,1:j);
            saved_iterations = saved_iterations(1:j);
            break
        end
    end
    if NE_Rnrm(k) <= NE_Rtol
        if verbose
            disp('Normal equations residual tolerance satisfied')
        end
        StopFlag = 'Normal equations residual tolerance satisfied';
        StopIt = k;
        StopReg.It = k;
        StopReg.X = x;
        if errornorms, StopReg.Enrm = Enrm(k); end
        if ~ NoStop
            if ~AlreadySaved
                j = j+1;
                X(:,j) = x;
                saved_iterations(j) = k;
                AlreadySaved = 1;
            end
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if errornorms, Enrm = Enrm(1:k); end
            X = X(:,1:j);
            saved_iterations = saved_iterations(1:j);
            break
        end 
    end
end
if k == MaxIter
  if StopIt == MaxIter
    % Stop because max number of iterations reached
    if verbose
        disp('Reached maximum number of iterations')
    end
    StopFlag = 'Reached maximum number of iterations';
    StopReg.It = k;
    StopReg.X = x;
    if errornorms, StopReg.Enrm = Enrm(k); end
    if ~AlreadySaved
        j = j+1;
        X(:,j) = x;
        saved_iterations(j) = k;
    end
    Xnrm    = Xnrm(1:k);
    Rnrm    = Rnrm(1:k);
    NE_Rnrm = NE_Rnrm(1:k);
    if errornorms, Enrm = Enrm(1:k); end
    X = X(:,1:j);
    saved_iterations = saved_iterations(1:j);
  end 
end
if ~noIterBar, close(h_wait), end
if nargout==2
  info.its = k;
  info.saved_iterations = saved_iterations(1:j);
  info.StopFlag = StopFlag;
  info.StopReg = StopReg;
  info.Rnrm = Rnrm(1:k);
  info.NE_Rnrm = NE_Rnrm(1:k);
  info.Xnrm = Xnrm(1:k);
  if errornorms
    info.Enrm = Enrm(1:k);
    info.BestReg = BestReg;
  end
end