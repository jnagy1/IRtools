function [X,info] = IRnnfcgls(A,b,varargin)
% IRnnfcgls Modified flexible CGLS for nonnegatively constrained LS problems
%
% options  = IRnnfcgls('defaults')
% [X,info] = IRnnfcgls(A,b)
% [X,info] = IRnnfcgls(A,b,K)
% [X,info] = IRnnfcgls(A,b,options)
% [X,info] = IRnnfcgls(A,b,K,options)
%
% This function implements a modified flexible CGLS method in which the
% preconditioner enforces nonnegativity.
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
%  K : (optional) integer vector that specifies which (total) iterates are 
%      returned in X; the maximum number of iterations is assumed to be max(K)
%      [ positive integer | vector of positive components ]
%  options : structure with the following fields (optional)
%      x0         - initial guess for the iterations; default = zero vector
%                   [ array | b | {'none'} ]
%      MaxIterIn  - maximum number of inner iterations
%      MaxIterOut - maximum number of outer iterations
%      x_true     - true solution; allows us to returns error norms with
%                   respect to x_true at each iteration
%                   [ array | {'none'} ]
%      stopCrit   - stopping criterion for the (total) iterations
%                   [ {'none'} | 'discrep' ]
%                   Note: 'discrepancy' requires NoiseLevel and eta
%      NoiseLevel - norm of noise in rhs divided by norm of rhs
%                   [ {'none'} | nonnegative scalar ]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      trunc      - truncation parameter for the FCGLS direction update
%                   [ positive integer | {'none'} | 'full' ]
%      tollalpha  - if the steplength of the CGLS update is less than
%                   tollalpha, no update happens and the inner
%                   iterations are terminated
%                   [ {1e-15} | positive scalar (close to 0) ] 
%      IterBar    - shows the progress of the iterations
%                   [ {'on'} | 'off' ]
%      NoStop     - specifies whether the iterations should proceed
%                   after a stopping criterion has been satisfied
%                   [ 'on' | {'off'} ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its      - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - string that describes the output/stopping condition:
%                   * Immediate stagnation of the method (in this case,
%                     a different x0 should be assigned)
%                   * Stagnation of the method
%                   * Reached maximum number of iterations
%                   * Discrepancy principle satisfied
%      StopReg  - struct containing information about the solution that
%                 satisfies the stopping criterion, with the fields:
%                   It   : iteration where the stopping criterion is satisfied
%                   X    : the solution satisfying the stopping criterion
%                   Enrm : the corresponding relative error (requires x_true)
%      Rnrm     - relative residual norms at each iteration
%      NE_Rnrm  - normal eqs relative residual norms at each iteration
%      Xnrm     - solution norms at each iteration
%      Enrm     - relative error norms (requires x_true) at each iteration
%      Xout     - approximate solutions at the beginning of each inner
%                 cycle (stored column-wise)
%      itsInOut - 3-column matrix, where the columns store
%                   1. outer iteration count
%                   2. inner iteration count (i.e., for each cycle)
%                   3. total iteration count
%      Steps    - steplength for the CGLS update at each iteration
%
% See also: IRcgls, IRconstr_ls, IRfista, IRmrnsd, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('trunc', 'none', 'x0', 'none', 'MaxIterIn', 10, ...
    'MaxIterOut', 20 , 'stopCrit', 'none', 'x_true', 'none', ...
    'IterBar', 'on', 'NoStop', 'off', 'verbosity', 'on', ...
    'NoiseLevel', 'none', 'eta', 1.01, 'tollalpha', 1e-15);
  
% If input is 'defaults,' return the default options in X.
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    X = defaultopt;
    return;
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

trunc      = IRget(options, 'trunc',     [], 'fast');
MaxIterIn  = IRget(options, 'MaxIterIn', [], 'fast');
MaxIterOut = IRget(options, 'MaxIterOut',[], 'fast');
x_true     = IRget(options, 'x_true',    [], 'fast');
stopCrit   = IRget(options, 'stopCrit',  [], 'fast');
IterBar    = IRget(options, 'IterBar',   [], 'fast');
NoStop     = IRget(options, 'NoStop',    [], 'fast');
NoiseLevel = IRget(options, 'NoiseLevel',[], 'fast');
eta        = IRget(options, 'eta',       [], 'fast');
tollalpha  = IRget(options, 'tollalpha', [], 'fast');
verbose    = IRget(options, 'verbosity', [], 'fast');

NoStop = strcmp(NoStop,'on');
verbose = strcmp(verbose, 'on');

if isscalar(NoiseLevel)
    stopCrit = 'discrep';
end

if strcmp(stopCrit, 'none')
    NoStop = 1;
end

MaxIter = MaxIterIn*MaxIterOut;
% Note: MaxIter can be quite big.  Anyway, if no stopping criterion applies
% or if NoStop = 'on', our choice is to perform MaxIter total iterations
% and warn the user about this

% Setting K.
if isempty(K)
    K = MaxIter;
    if NoStop
        warning(['Note that with options.NoStop = ''on'' or without any specified stopping criterion ',...
                 'the maximum number of ', num2str(MaxIter), ' iterations are performed; ', ...
                 'to reduce this number specify/change K or options.MaxIterIn and options.MaxIterOut'])
    end
end

% Checking the entries in K, and sorting the iterations (in case they are
% shuffled in input).
K = K(:); K = sort(K,'ascend'); K = unique(K);
if ~((isreal(K) && (all(K > 0)) && all(K == floor(K))))
    error('K must be a vector of positive real integers')
end
if K(end) ~= MaxIter
    MaxIter = K(end);    
end

if ischar(NoiseLevel) && strcmp(stopCrit,'discrep')
    error('The noise level must be assigned')
end

if (ischar(trunc) && strcmp(trunc, 'none'))
    trunc = MaxIterIn;
elseif (ischar(trunc) && strcmp(trunc, 'full'))
    trunc = 0;
end

z = Atransp_times_vec(A, b);
n = length(z(:));
nrmb = norm(b(:));
nrmz = norm(z(:));

% Compute the product A*x0, which is useful for setting the initial residual.
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

x(x<0) = 0;
Ax = A_times_vec(A, x); Ax = Ax(:);
r  = b(:) - Ax(:);
z = Atransp_times_vec(A, r);

% Declare matrices.
X       = zeros(n,length(K));
Xout    = zeros(n,MaxIterOut);
Xnrm    = zeros(MaxIter,1);
Rnrm    = zeros(MaxIter,1);
NE_Rnrm = zeros(MaxIter,1);
IterIn  = zeros(MaxIterOut,3);
Alpha   = zeros(MaxIter,1);
saved_iterations = zeros(1, length(K));
StopReg.It = MaxIter;
StopReg.X = [];
if strcmp(x_true,'none')
    errornorms = false;
else
    errornorms = true;
    nrmx = norm(x_true(:));
    Enrm = zeros(max(K),1);
    BestReg.It = [];
    BestReg.X = [];
    BestReg.Enrm = [];
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestReg.NE_Rnrm = [];
    BestEnrm = 1e10;
end

j   = 0; % Counter for the stored iterations (in X and saved_iterations)
k   = 0; % Counter for the total iterations
kout= 0; % Counter for the outer iterations
stagn = 0; % Detects if the method stagnates

% param is a struct containing the inputs and the outputs of the inner
% iterations (i.e., it is updated within the outer cycles).
param.MaxIter   = MaxIter;
param.MaxIterIn = MaxIterIn;
param.trunc = trunc;
param.tollalpha = tollalpha;
param.K = K;
param.k = k;
param.j = j;
param.r = r;
param.z = z(:);
param.X = X;
param.Xnrm = Xnrm;
param.Rnrm = Rnrm;
param.NE_Rnrm = NE_Rnrm;
param.Steps = Alpha;
param.errornorms = errornorms;
param.x_true = x_true(:);
if errornorms
    param.Enrm = Enrm;
    param.BestEnrm = BestEnrm;
    param.BestReg = BestReg;
    param.nrmx = nrmx;
end
param.saved_iterations = saved_iterations;
param.exitflag = 'none';
param.Ax = Ax(:);
param.NoStop = NoStop;
param.stopCrit = stopCrit;
param.StopReg = StopReg;
param.eta = eta;
param.NoiseLevel = NoiseLevel;
param.stagn = stagn;
param.verbose = verbose;
param.nrmb = nrmb;
param.nrmz = nrmz;

% Iterate.
noIterBar = strcmp(IterBar, {'off'});
if ~noIterBar
  h_wait = waitbar(0, 'Running iterations, please wait ...');
end

% Potentially iterate till the maximum number of (total) iterations
% is reached.
while k < MaxIter 
    kout = kout + 1;
    x(x<0) = 0; % Redundant, just to make sure nonnegativity is respected.
    Xout(:,kout) = x(:);
    if ~noIterBar
        waitbar(k/MaxIter, h_wait)
    end
    [x, iterin, k, param] = nnfcgls_in(A, b, x, k, param);
    IterIn(kout,:) = [kout, iterin];
    if (~ (NoStop) && ~ (strcmp(param.exitflag, 'none'))) || (param.stagn)
    % This means that:
    % 1. a stopping rule  as been reached, and we indeed want to stop; OR
    % 2. the method stagnates.
    % Therefore we also break the outer iterations.
        break
    end  
end
if k == MaxIter && strcmp(param.exitflag, 'none')
    param.exitflag = 2;
    if verbose
        disp('Reached maximum number of iterations')
    end
    param.StopReg.It = k;
    param.StopReg.X = x(:,end);
    if errornorms, param.StopReg.Enrm = param.Enrm(k); end
end
switch param.exitflag
    case 0
        StopFlag = 'Immediate stagnation of the method';
    case 1
        StopFlag = 'Stagnation of the method';
    case 2
        StopFlag = 'Reached maximum number of (total) iterations';
    case 3
        StopFlag = 'Discrepancy principle satisfied';
end
IterIn = IterIn(1:kout, :);
if ~noIterBar
    waitbar(k/MaxIter, h_wait);
    close(h_wait)
end
j    = param.j;
X    = param.X;
X    = X(:,1:j);
if isempty(X)
    X = x;
end
Xout = Xout(:,1:kout);
if nargout==2
  info.its = k;
  info.StopReg = param.StopReg;
  info.Xout = Xout;
  info.itsInOut = IterIn;
  info.saved_iterations = param.saved_iterations;
  info.Steps = param.Steps;
  info.Rnrm = param.Rnrm;
  info.NE_Rnrm = param.NE_Rnrm;
  info.Xnrm = param.Xnrm;
  if errornorms
    info.Enrm = param.Enrm;
    info.BestReg = param.BestReg;
  end
  info.StopFlag = StopFlag;
end

% AUXILIARY FUNCTIONS ------------------------------------------------

function [x, iterin, k, param] = nnfcgls_in(A, b, x, k, param)
% These parameters are not updated within the cycles.
MaxIter = param.MaxIter;
MaxIterIn = param.MaxIterIn;
trunc = param.trunc;
K = param.K;
tollalpha = param.tollalpha;
errornorms = param.errornorms;
x_true = param.x_true;
if errornorms, nrmx = param.nrmx; end
nrmb = param.nrmb;
nrmz = param.nrmz;
NoStop = param.NoStop;
stopCrit = param.stopCrit;
eta = param.eta;
NoiseLevel = param.NoiseLevel;
verbose = param.verbose;
% These parameters/variables need to be updated.
X = param.X;
j = param.j;
r = param.r;
z = param.z;
Alpha = param.Steps;
exitflag = param.exitflag;
saved_iterations = param.saved_iterations;
Xnrm = param.Xnrm;
Rnrm = param.Rnrm;
NE_Rnrm = param.NE_Rnrm;
if errornorms
    Enrm = param.Enrm;
    BestEnrm = param.BestEnrm;
    BestReg = param.BestReg;
end
StopReg = param.StopReg;
StopIt = param.StopReg.It;
stagn = param.stagn;
if trunc ~= 0
    Q = zeros(size(x,1), MaxIterIn+1);
    AQ = zeros(size(r,1), MaxIterIn+1);
end
%
if max(x(:))~=0
    zb = x(:).*z(:);
else
    zb = z(:);
end
q = zb;
Azb = A_times_vec(A, zb);
Aq = Azb;
if trunc
    Q(:,1) = q(:);
    AQ(:,1) = Aq(:);
end
% Iterate
for i = 1:MaxIterIn
    k = k+1; 
    w = Aq;
    theta = (r(:)'*w(:))/((w(:))'*w(:));
    neg_ind = q < 0;
    alpha = min( theta, min( -x(neg_ind)./q(neg_ind) ) );
    % Checks.
    if isempty(alpha)
        alpha = theta;
    end
    if (alpha < tollalpha) && (i>1)
        iterin = [i-1, k-1];
        k = k-1;
        break
    elseif (alpha < tollalpha) && (i==1) && (k==1)
        disp('The method immediately stagnates; try a different initial guess')
        stagn = 1;
        iterin = [0, 0];
        StopReg.It = 0;
        StopReg.X = [];
        exitflag = 0;
        Alpha = [];
        Rnrm = [];
        NE_Rnrm = [];
        Xnrm = [];
        if errornorms
            Enrm = [];
            StopReg.Enrm = [];
        end
        break
    elseif (alpha < tollalpha) && (i==1) && (k>1)
        disp('The method stagnates; no further improvements are possible')
        stagn = 1;
        iterin = [0, k-1];
        k = k-1;
        if StopIt == MaxIter
            StopReg.It = k;
            StopReg.X = x(:);
            exitflag = 1;
        end
        Alpha = Alpha(1:k);
        Rnrm = Rnrm(1:k);
        NE_Rnrm = NE_Rnrm(1:k);
        Xnrm = Xnrm(1:k);
        if errornorms
            Enrm = Enrm(1:k);
            StopReg.Enrm = Enrm(k);
        end
        break
    end
    if k > MaxIter
        if StopIt == MaxIter
            if verbose
                disp('Reached maximum number of iterations')
            end
            StopReg.It = k-1;
            StopReg.X = x(:);
            if errornorms, StopReg.Enrm = Enrm(k-1); end
            exitflag = 2;
        end
        iterin = [i-1, k-1];
        k = k-1;
        break
    end
    % Compute the next approximation.
    Alpha(k) = alpha;
    AlreadySaved = 0;
    x = x + alpha*q;
    Xnrm(k) = norm(x(:));
    r = r - alpha*w;
    Rnrm(k) = norm(r(:))/nrmb;
    z = Atransp_times_vec(A, r);
    NE_Rnrm(k) = norm(z(:))/nrmz;
    zb = x(:).*z(:); 
    Azb = A_times_vec(A, zb);
    if trunc > 0
        beta = zeros(i,1);
        for l = 1:i
            beta(l) = -(Azb(:)'*AQ(:,l))/(AQ(:,l)'*AQ(:,l));
        end
        if trunc >= i
                q = zb + Q(:,1:i)*beta(1:i);
                Aq = Azb + AQ(:,1:i)*beta(1:i);
        else
            trunctemp = i-trunc;
            q = zb + Q(:,trunctemp+1:trunctemp+trunc)*beta(trunctemp+1:trunctemp+trunc);
            Aq = Azb + AQ(:,trunctemp+1:trunctemp+trunc)*beta(trunctemp+1:trunctemp+trunc);
        end
        Q(:,i+1) = q(:);
        AQ(:,i+1) = Aq(:);
    else
        beta = -(Azb(:)'*w(:))/(w(:)'*w(:));
        q = zb + beta*q;
        Aq = Azb + beta*Aq;
    end
    
    if any(k==K) && ~ AlreadySaved
        j = j+1;
        saved_iterations(j) = k;
        X(:,j) = x;
        AlreadySaved = 1;
    end
    
    if errornorms
        Enrm(k) = norm(x(:) - x_true(:))/nrmx;
        if Enrm(k)<BestEnrm
            BestReg.It = k;
            BestReg.X = x(:);
            BestEnrm = Enrm(k);
            BestReg.Enrm = BestEnrm;
            BestReg.Xnrm = Xnrm(k);
            BestReg.Rnrm = Rnrm(k);
            BestReg.NE_Rnrm = NE_Rnrm(k);
        end
    end
    
    if strcmp(stopCrit, 'discrep')
        if Rnrm(k) < eta*NoiseLevel
            if StopIt >= k
                if verbose
                    disp('The discrepancy principle is satisfied')
                end
                exitflag = 3;
                StopReg.It = k;
                StopReg.X = x;
                if errornorms, StopReg.Enrm = Enrm(k); end
                StopIt = k;
            end
            if ~ NoStop
                Xnrm    = Xnrm(1:k);
                Rnrm    = Rnrm(1:k);
                NE_Rnrm = NE_Rnrm(1:k);
                Alpha = Alpha(1:k);
                if errornorms, Enrm = Enrm(1:k); end
                if ~AlreadySaved
                    j = j+1;
                    saved_iterations(j) = k;
                    X(:,j) = x;
                end
                X = X(:,1:j);
                saved_iterations = saved_iterations(1:j);
                iterin = [i, k];
                break
            end
        end
    end
    if i == MaxIterIn
        iterin = [i, k];
    end
end
param.X = X;
param.j = j;
param.r = r;
param.z = z;
param.Steps = Alpha;
param.saved_iterations = saved_iterations;
param.Xnrm = Xnrm;
param.Rnrm = Rnrm;
param.NE_Rnrm = NE_Rnrm;
if errornorms
    param.Enrm = Enrm;
    param.BestEnrm = BestEnrm;
    param.BestReg = BestReg;
end
param.Ax = b(:) - r(:);
param.StopReg = StopReg;
param.stagn = stagn;
param.exitflag = exitflag;