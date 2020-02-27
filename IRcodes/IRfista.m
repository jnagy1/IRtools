function [X,info] = IRfista(A,b,varargin)
%IRfista FISTA algorithm for constrained least squares problems
%
% options  = IRfista('defaults')
% [X,info] = IRfista(A,b)
% [X,info] = IRfista(A,b,K)
% [X,info] = IRfista(A,b,options)
% [X,info] = IRfista(A,b,K,options)
%
% This function uses the first-order optimization method FISTA to solve
% the solves the least squares or Tikhonov problem with constraints.
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
%      x0         - initial guess for the iterations; default = zero vector
%                   [ array | {'none'} ]
%      MaxIter    - maximum allowed number of CGLS iterations
%                   [ positive integer | {100}]
%                   NOTE: K overrules MaxIter if both are assigned
%      x_true     - true solution; allows us to returns error norms with
%                   respect to x_true at each iteration
%                   [ array | {'none'} ]
%      NoiseLevel - norm of noise in rhs divided by norm of rhs 
%                   [ {'none'} | nonnegative scalar]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      NE_Rtol    - relative tolerance on the normal eqs residual norm
%                   [ {1e-12} | positive integer ]
%      NoStop     - specifies whether the iterations should proceed after
%                   a stopping criterion has been satisfied
%                   [ 'on' | {'off'}]
%      RegParam   - regularization parameter lambda, to be employed 
%                   if FISTA is used to solve the regularized problem
%                   (A'*A + lambda^2*L'*L)*x = A'*b;
%                   [ nonnegative scalar | {0} ]
%      xMin       - lower bound for the solution elements
%                   [ {0} | scalar value ]
%      xMax       - upper bound for the solution elements
%                   [ {Inf} | scalar value ]
%      xEnergy    - value of the energy constraint
%                   [ {'none'} | positive scalar ]
%      shrink     - determine wether iterative shrinkage thresholding is applied
%                   [ {'off'} | 'on' ]
%      IterBar    - shows the progress of the iterations
%                   [ {'on'} | 'off' ]
%   SparsityTrans - sparsity transform for the solution
%                   [ {'none'} | 'dwt' ]
%      wname      - discrete wavelet transform name (meaningful if 
%                   SpartistyTrans is 'dwt')
%                   [ {'db1'} ]
%      wlevels    - discrete wavelet transform level (meaningful if 
%                   SpartistyTrans is 'dwt')
%                   [ {2} | positive integer]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its      - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - a string that describes the stopping condition:
%                   * Performed max number of iterations
%                   * Residual tolerance satisfied (discrepancy principle) 
%                   * Normal equationa residual tolerance satisfied
%      Rnrm     - relative residual norms at each iteration
%      NE_Rnrm  - normal eqs relative residual norms at each iteration
%      Xnrm     - solution norms at each iteration
%      Enrm     - relative error norms (requires x_true) at each iteration
%      StopReg  - struct containing information about the solution that
%                 satisfies the stopping criterion.  Fields:
%                   It :  iteration where the stopping criterion is satisfied
%                   X  :  solution satisfying the stopping criterion
%                   Enrm: the corresponding relative error (requires x_true)
%      BestReg  - struct containing information about the solution that
%                 minimizes Enrm (requires x_true).  Fields:
%                   It :  iteration where the minimum is attained
%                   X  :  best solution
%                   Enrm: best relative error
%
% See also: IRconstr_ls, IRmrnsd, IRnnfcgls, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0','none', 'MaxIter',100 , 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, 'IterBar','on', ...
    'NoStop','off', 'shrink', 'off', 'RegParam',0, 'xMin',0, 'xMax',Inf,...
    'xEnergy','none','stepsize','none',...
    'backtracking','off','backit',10,'backscalar',1.1,...
    'SparsityTrans', 'none', 'wname', 'db1', 'wlevels', 2);
  
% If input is 'defaults,' return the default options in X
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

MaxIter    = IRget(options, 'MaxIter',      [], 'fast');
x_true     = IRget(options, 'x_true',       [], 'fast');
NoiseLevel = IRget(options, 'NoiseLevel',   [], 'fast');
eta        = IRget(options, 'eta',          [], 'fast');
NE_Rtol    = IRget(options, 'NE_Rtol',      [], 'fast');
IterBar    = IRget(options, 'IterBar',      [], 'fast');
TikParam   = IRget(options, 'RegParam',     [], 'fast');
xMin       = IRget(options, 'xMin',         [], 'fast');
xMax       = IRget(options, 'xMax',         [], 'fast');
xEnergy    = IRget(options, 'xEnergy',      [], 'fast');
NoStop     = IRget(options, 'NoStop',       [], 'fast');
shrinkage  = IRget(options, 'shrink',       [], 'fast');
t          = IRget(options, 'stepsize',     [], 'fast');
backtrack  = IRget(options, 'backtracking', [], 'fast');
SparsTrans = IRget(options, 'SparsityTrans',[], 'fast');


NoStop = strcmp(NoStop,'on'); 

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

if ischar(t)
    normestA = IRnormest(A, b);
    t = 1/(normestA^2 + TikParam^2);
end

if isempty(NoiseLevel) || strcmp(NoiseLevel,'none')
    Rtol = 0;
else
    Rtol = eta*NoiseLevel;
end

shrinkage = strcmp(shrinkage, 'on');
if shrinkage
    if ischar(TikParam) || (isscalar(TikParam) && TikParam == 0)
        error('A positive value for the regularization parameter should be assigned')
    else
        thrshrink = TikParam*t;
    end
end

backtrack = strcmp(backtrack, 'on');
if backtrack
    backit     = IRget(options, 'backit', [], 'fast');
    backscalar = IRget(options, 'backscalar', [], 'fast');
    if backscalar < 1
        error(' The backtracking constant should be strictly greater than 1')
    end
end

%  We need to find the number of columns in matrix A, but if A is not given 
%  as a matrix, and no initial guess is given, then we can find it by 
%  computing A'*b.  We need this anyway, so doesn't cost additional work.
d = Atransp_times_vec(A, b);
n = length(d);

if strcmp(SparsTrans, 'none')
    Trans = speye(n);
elseif strcmp(SparsTrans, 'dwt')
    wname   = IRget(options, 'wname',   [], 'fast');
    wlevels = IRget(options, 'wlevels', [], 'fast'); 
    if wlevels > 1/2*(log2(n))
        error('The assigned wavelet levels are too high. Make sure that wlevels <= 1/2*(log2(n))')
    end
    Trans = FreqMatrix('dwt', [sqrt(n) sqrt(n)], wname, wlevels);
end

% See if an initial guess is given; if not use 0 as the initial guess.  
x = IRget(options, 'x0', [], 'fast');

if strcmp(x,'none')
    if ~shrinkage
        % the default initial guess for the iterations is defined as the
        % minimizer of || b - A*(alpha*ones(n,1)) ||_2
        coeffx0 = A_times_vec(A, ones(n,1)); alpha = (coeffx0'*b)/norm(coeffx0)^2;
        if coeffx0 <= 0
            alpha = sqrt(eps);
        end
        x = alpha*ones(n,1);
    else
        x = zeros(n,1);
    end
end

Ax = A_times_vec(A, x);
r = b - Ax;
d = Atransp_times_vec(A, r);

% Declare matrices.
X = zeros(n,length(K));
saved_iterations = zeros(1, length(K));
Xnrm    = zeros(MaxIter,1);
Rnrm    = zeros(MaxIter,1);
NE_Rnrm = zeros(MaxIter,1);
if strcmp(x_true,'none')
    errornorms = false;
else
    errornorms = true;
    Enrm = zeros(max(K),1);
    nrmtrue = norm(x_true(:));
    BestReg.It = [];
    BestReg.X =[];
    BestReg.Enrm = [];
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestReg.NE_Rnrm = [];
    BestEnrm = 1e10;
end
nrmb = norm(b(:));
nrmAtb = norm(d(:));

%  We need to initialize this method with the first two iterations.

%  Here is iteration 1:
if shrinkage 
    xnew = x + t*d;
    try
        xnew = Trans*xnew;
        xnew = Shrink(xnew, thrshrink);
        xnew = Trans'*xnew;
    catch
        error('Check the length of x0')
    end    
else
    xnew = x + t*(d - TikParam^2*x);
end
xnew = Project(xnew, xMin, xMax, xEnergy);

% stuff needed for backtracking
% Ax = A_times_vec(A, x);
% r = b - Ax;
% 
Axnew = A_times_vec(A, xnew);
rnew = b - Axnew;
if backtrack
    if shrinkage
        Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew,1);
        Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew,1);
    else
        Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew)^2;
        Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew)^2;
    end
    backittemp = 1;
    while Fx>Qxy && backittemp < backit
        t = t/backscalar;
        if shrinkage
            thrshrink = TikParam*t;
            xnew = x + t*d;
            xnew = Trans*xnew; xnew = Shrink(xnew, thrshrink); xnew = Trans'*xnew;
        else
            xnew = x + t*(d - TikParam^2*x);
        end
        xnew = Project(xnew, xMin, xMax, xEnergy);
        Axnew = A_times_vec(A, xnew);
        rnew = b - Axnew;
        if shrinkage
            Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew,1);
            Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew,1);
        else
            Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew)^2;
            Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew)^2;
        end
        backittemp = backittemp + 1;
    end
end
x = xnew;
r = rnew;
x_save1 = x;
j = 0;    
if any(K == 1)
    j = j+1;
    X(:,j) = x;
    saved_iterations(j) = 1;
end

Rnrm(1) = norm(r)/nrmb;
d = Atransp_times_vec(A, r);
NE_Rnrm(1) = norm(d)/nrmAtb;
if errornorms
    Enrm(1) = norm(x_true-x)/nrmtrue;
    if Enrm(1)<BestEnrm
        BestReg.It = 1;
        BestReg.X = x;
        BestEnrm = Enrm(1);
        BestReg.Enrm = BestEnrm;
        BestReg.Xnrm = Xnrm(1);
        BestReg.Rnrm = Rnrm(1);
        BestReg.NE_Rnrm = NE_Rnrm(1);
    end
end

% And here is iteration 2:
tk = 0.5*(1+sqrt(5));
if shrinkage 
    xnew = x + t*d;
    xnew = Trans*xnew;
    xnew = Shrink(xnew, thrshrink);
    xnew = Trans'*xnew;
else
    xnew = x + t*(d - TikParam^2*x);
end
xnew = Project(xnew, xMin, xMax, xEnergy);
Axnew = A_times_vec(A, xnew);
rnew = b - Axnew;
if backtrack
    if shrinkage
        Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew,1);
        Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew,1);
    else
        Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew)^2;
        Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew)^2;
    end
    backittemp = 1;
    while Fx>Qxy && backittemp < backit
        t = t/backscalar;
        if shrinkage
            thrshrink = TikParam*t;
            xnew = x + t*d;
            xnew = Trans*xnew; xnew = Shrink(xnew, thrshrink); xnew = Trans'*xnew;
        else
            xnew = x + t*(d - TikParam^2*x);
        end
        xnew = Project(xnew, xMin, xMax, xEnergy);
        Axnew = A_times_vec(A, xnew);
        rnew = b - Axnew;
        if shrinkage
            Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew,1);
            Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew,1);
        else
            Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew)^2;
            Qxy = 1/2*norm(r)^2 + (xnew - x)'*d(:) + 1/(2*t)*norm(x - xnew)^2 + TikParam*norm(xnew)^2;
        end
        backittemp = backittemp + 1;
    end
end
x = xnew;
r = rnew;
x_save2 = x;
if any(K == 2)
    j = j+1;
    X(:,j) = x;
    saved_iterations(j) = 2;
end
Rnrm(2) = norm(r)/nrmb;
d = Atransp_times_vec(A, r);
NE_Rnrm(2) = norm(d)/nrmAtb;
if errornorms
    Enrm(2) = norm(x_true-x)/nrmtrue;
    if Enrm(2)<BestEnrm
        BestReg.It = 2;
        BestReg.X = x;
        BestEnrm = Enrm(2);
        BestReg.Enrm = BestEnrm;
        BestReg.Xnrm = Xnrm(2);
        BestReg.Rnrm = Rnrm(2);
        BestReg.NE_Rnrm = NE_Rnrm(2);
    end
end

% Iterate.
noIterBar = strcmp(IterBar,{'off'});
if ~noIterBar
  h_wait = waitbar(0, 'Running iterations, please wait ...');
end
for k=3:MaxIter
    AlreadySaved = 0;
    if ~noIterBar
        waitbar(k/MaxIter, h_wait)
    end
    tk_save = tk;
    tk = 0.5*(1 + sqrt(1 + 4*tk_save^2));
    y = x_save2 + ((tk_save-1)/tk)*(x_save2 - x_save1);
    Ay = A_times_vec(A, y);
    d = Atransp_times_vec(A, b-Ay);
    if shrinkage 
        xnew = y + t*d;
        xnew = Trans*xnew;
        xnew = Shrink(xnew, thrshrink);
        xnew = Trans'*xnew;
    else
        xnew = y + t*(d - TikParam^2*y);
    end
    xnew = Project(xnew, xMin, xMax, xEnergy);
    Axnew = A_times_vec(A, xnew);
    rnew = b - Axnew;
    if backtrack
        if shrinkage
            Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew,1);
            Qxy = 1/2*norm(r)^2 + (xnew - y)'*d(:) + 1/(2*t)*norm(xnew - y)^2 + TikParam*norm(xnew,1);
        else
            Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew)^2;
            Qxy = 1/2*norm(r)^2 + (xnew - y)'*d(:) + 1/(2*t)*norm(xnew - y)^2 + TikParam*norm(xnew)^2;
        end
        backittemp = 1;
        while Fx>Qxy && backittemp < backit
            t = t/backscalar;
            if shrinkage
                thrshrink = TikParam*t;
                xnew = y + t*d;
                xnew = Trans*xnew; xnew = Shrink(xnew, thrshrink); xnew = Trans'*xnew;
            else
                xnew = y + t*(d - TikParam^2*y);
            end
            xnew = Project(xnew, xMin, xMax, xEnergy);
            Axnew = A_times_vec(A, xnew);
            rnew = b - Axnew;
            if shrinkage
                Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew,1);
                Qxy = 1/2*norm(r)^2 + (xnew - y)'*d(:) + 1/(2*t)*norm(xnew - y)^2 + TikParam*norm(xnew,1);
            else
                Fx = 1/2*norm(rnew)^2 + TikParam*norm(xnew)^2;
                Qxy = 1/2*norm(r)^2 + (xnew - y)'*d(:) + 1/(2*t)*norm(xnew - y)^2 + TikParam*norm(xnew)^2;
            end
            backittemp = backittemp + 1;
        end
    end
    x = xnew;
    r = rnew;
    if any(k==K) && ~AlreadySaved
        j = j+1;
        X(:,j) = x;
        saved_iterations(j) = k;
        AlreadySaved = 1; 
    end
    x_save1 = x_save2;
    x_save2 = x;
    % Compute norms.
    Xnrm(k) = norm(x);
    Rnrm(k) = norm(r)/nrmb;
    NE_Rnrm(k) = norm(Atransp_times_vec(A,r))/nrmAtb;
    
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
    if Rnrm(k) <= Rtol  && (StopIt == MaxIter)
        disp('Residual tolerance satisfied')
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
    if NE_Rnrm(k) <= NE_Rtol && (StopIt == MaxIter)
        disp('Normal equations residual tolerance satisfied')
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
        % Stop because max number of iterations reached.
        disp('Reached maximum number of iterations')
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
  info.Rnrm = Rnrm(1:k);
  info.NE_Rnrm = NE_Rnrm(1:k);
  info.Xnrm = Xnrm(1:k);
  if errornorms
    info.Enrm = Enrm(1:k);
    info.BestReg = BestReg;
  end
  info.StopFlag = StopFlag;
  info.StopReg = StopReg;
  info.saved_iterations = saved_iterations(1:j);
end


function x = Project(x, xMin, xMax, xEnergy)
%
%  Compute the projection
%
if strcmpi(xEnergy,'none')
    x = min(x,xMax);
    x = max(x,xMin);
else
    x = gdnnf_projection(x, xEnergy);
end

function x = Shrink(x, T)
%
%  Compute the shrinkage
%
x = sign(x).*max(abs(x)-T, 0);
