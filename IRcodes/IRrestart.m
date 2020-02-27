function [X, info] = IRrestart(A, b, varargin)
% IRrestart Restarted Krylov subspace methods
%
% options  = IRrestart('defaults')
% [X,info] = IRrestart(A,b)
% [X,info] = IRrestart(A,b,K)
% [X,info] = IRrestart(A,b,options)
% [X,info] = IRrestart(A,b,K,options)
%
% This function provides a general framework for a variety of methods
% (IRconstr_ls, IRirn, IRhtv), and implements an inner-outer (or restarted)
% iteration scheme.
%
% Semi-convergent or hybrid iterative solvers are used in the inner iterations, 
% using one of the iterative methods in IRtools (e.g., IRhybrid_gmres).
% Every outer iteration produces a new approximation that incorporates the 
% desired properties and/or constraints (see IRconstr_ls, IRirn, IRhtv
% for further details).
%
% The regularization parameter and number of inner iterations influence
% the behavior and convergence of the outer iterations.
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
%  b : rhs vector
%  K : (optional) integer vector that specifies which (total) iterates are 
%      returned in X; the maximum number of iterations is assumed to be max(K)
%      [ positive integer | vector of positive components ]
%  options : structure with the following fields (optional)
%             x0 - initial guess for the iterations; default = zero vector
%                        [ array | {'none'} ]
%             MaxIterIn - maximum number of inner iterations
%             MaxIterOut - maximum number of outer iterations
%             x_true - true solution; allows us to returns error norms
%                      with respect to x_true at each iteration
%                      [ array | {'none'} ]
%             RegParam - a value or a method to find the regularization
%                        used in the inner iterations: 
%                        [non-negative scalar | {'gcv'} | 'discrep' ]
%                    This also determines which stopping rule is used for
%                    the inner iterations.
%                    If 'gcv' is chosen, the inner iteration is stopped when
%                    the GCV function minimum stabilizes or increases 
%                    within a certain window of iterations (see 'stopGCV',
%                    'FlatTol' and 'MinTol').
%                    If 'discrep' is chosen, and NoiseLevel is
%                    provided, then the discrepancy principle is used
%                    for a stopping criterion (see 'NoiseLevel' and 'eta').
%             stopGCV - stopping criterion for the inner iterations when
%                       GCV is used
%                       [ GCVvalues | {'resflat'} ]
%             FlatTol - tolerance for detecting flatness (stabilization)
%                       in the GCV function as a stopping criterion for the
%                       inner iterations
%                       [ non-negative scalar | {10^-6} ]
%             MinTol - window of iterations - if the GCV minimum continues
%                      to increase over this window, then the inner
%                      iterations are stopped:
%                      [ positive integer | {3}]
%             RegMatrix - priorconditioner for the inner iterations
%                         [ {'identity'} | square nonsingular matrix |
%                           function handle ]
%             NoiseLevel - norm of noise in rhs divided by norm of rhs
%                          (must be assigned in RegParam is 'discrep')
%                          [ {none} | nonnegative scalar ]
%             eta - safety factor for the discrepancy principle
%                   [ 1.01 | scalar greater than (and close to) 1 ]
%             RegParam0 - first regularization parameter, used only on the 
%                         very first iteration (needed if RegParam is 'discrep')
%                         [ {1} | positive scalar ]
%             stopOut - stopping criterion for the outer iterations;
%                       [ {'xstab'} | 'Lxstab' | 'regPstab' ]
%             inSolver - solver to be employed during the inner iterations
%                       [ {'gmres'} | 'lsqr' | 'fgmres' | 'rrgmres' | 'cgls' ]
%             adaptConstr - approximate constraint or regularization
%                           to be incorporated
%                           [ {'tv'} | 'nn' | 'tvnn' | 'sp' | 'spnn' |
%                             'box' | 'energy' | 'project' | 'none']
%               'box' requires upper and lower scalar bounds 'xMin' and 'xMax'
%               'energy' requires the positive scalar 'xEnergy'
%               'project' requires 'xMin', 'xMax' and 'xEnergy'
%             xMin - lower bound for the solution elements
%                    [ {0} | scalar value ]
%             xMax - upper bound for the solution elements
%                    [ {Inf} | scalar value]
%             xEnergy - value of the energy constraint
%                       [ {'none'} | positive scalar ]
%             IterBar - shows the progress of the outer iterations
%                       [ {'on'} | 'off' ]
%             NoStopIn - Specifies whether the inner iterations should proceed
%                        after a stopping criterion has been satisfied
%                        [ 'on' | {'off'} ]
%             NoStopOut - specifies whether the outer iterations should
%                         proceed after a stopping criterion has been satisfied
%                         [ 'on' | {'off'} ]
%             warmrestart - specifies wether an available approximation of
%                           the solution should be used as an initial guess.
%                           [ 'on' | {'off'} ]
%             verbosity - switch on or off the "verbosity" of the function
%                       [ {'on'} | 'off' ]
% Note: the options structure can be created using the function 'IRset'.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%          its - number of the last computed iteration
%          saved_iterations - iteration numbers of iterates stored in X 
%          StopFlag_in - string that describes the inner stopping condition:
%              * Stopping criterion of the inner iterations is never satisfied
%              * Stopping criterion is satisfied at least once during the
%                inner iterations
%          StopFlag_out - string that describes the outer stopping condition;
%                         depending on the inputs it can be one of the following:
%              * Outer stopping criterion is never satisfied
%              * Diagonal weighting matrix is numerically zero
%              * Solution stabilizes
%              * Transformed solution stabilizes
%              * Regularization parameter stabilizes
%          Rnrm - relative residual norms at each iteration
%          Xnrm - solution norms at each iteration
%          Enrm - relative error norms (requires x_true) at each iteration
%          StopReg - struct containing information about the solution that
%                    satisfies the stopping criterion.  Fields:
%                      It : iteration where stopping criterion is satisfied
%                      X : solution satisfying the stopping criterion
%                      Enrm : the corresponding relative error (requires x_true)
%          BestReg - struct containing information about the solution that
%                    minimizes Enrm (requires x_true). Fields:
%                      It : iteration where the minimum is attained
%                      X : best solution
%                      Enrm : best relative error
%          Xout - approximate solutions at the end of each inner cycle,
%                 stored column-wise
%          itsInOut - 3-column matrix whose the columns store
%                       1. outer iteration count
%                       2. inner iteration count (i.e., for each cycle)
%                       3. total iteration count
%
% See also: IRconstr_ls, IRhtv, IRirn, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% September, 2017.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0', 'none', 'MaxIterIn', 30 , 'MaxIterOut', 20 , ...
    'RegParam', 'gcv', 'RegMatrix', 'identity', 'stopGCV', 'GCVvalues', ...
    'resflatTol', 0.05, 'GCVflatTol', 1e-6, 'GCVminTol', 3,...
    'x_true', 'none', 'IterBar', 'on', 'NoStop', 'off', 'NoStopIn', 'off', ...
    'NoStopOut', 'off', 'stopOut', 'xstab', 'stabOut', 1e-6, 'thr0', 1e-10, ...
    'NoiseLevel', 'none', 'eta', 1.01, 'RegParam0', 1,...
    'inSolver', 'gmres', 'adaptConstr', 'tv', 'verbosity', 'off',...
    'SparsityTrans', 'none', 'wname', 'db1', 'wlevels', 2, 'warmrestart', 'on');
  
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

MaxIterIn  = IRget(options, 'MaxIterIn',  [], 'fast');
MaxIterOut = IRget(options, 'MaxIterOut', [], 'fast');
RegParam   = IRget(options, 'RegParam',   [], 'fast');
L          = IRget(options, 'RegMatrix',  [], 'fast');
NoiseLevel = IRget(options, 'NoiseLevel', [], 'fast');
eta        = IRget(options, 'eta',        [], 'fast');
RegParamk  = IRget(options, 'RegParam0',  [], 'fast');
thr0       = IRget(options, 'thr0',       [], 'fast');
x_true     = IRget(options, 'x_true',     [], 'fast');
IterBar    = IRget(options, 'IterBar',    [], 'fast');
NoStop     = IRget(options, 'NoStop',     [], 'fast');
NoStopIn   = IRget(options, 'NoStopIn',   [], 'fast');
NoStopOut  = IRget(options, 'NoStopOut',  [], 'fast');
stopOut    = IRget(options, 'stopOut',    [], 'fast');
stabOut    = IRget(options, 'stabOut',    [], 'fast');
inSolver   = IRget(options, 'inSolver',   [], 'fast');
adaptConstr= IRget(options, 'adaptConstr',[], 'fast');
verbose    = IRget(options, 'verbosity',  [], 'fast');
warmrestart= IRget(options, 'warmrestart',[], 'fast');
SparsTrans = IRget(options, 'SparsityTrans', [], 'fast');

hybrid = (strcmp(inSolver, 'gmres') || strcmp(inSolver, 'fgmres') || strcmp(inSolver, 'lsqr'))...
    && (strcmp(adaptConstr, 'nn'));

NoStopIn = strcmp(NoStopIn,'on');
if strcmp(NoStop,'on')
    if ~NoStopIn
        warning(['The inner iterations are stopped when the stopping criterion is satisfied; ',...
            'to avoid this set ''NoStopIn'' to ''on'''])
        NoStop = 'off';
    else
        NoStop = 'on';
    end
end
if NoStopIn
    NoStop = 'on';
end
NoStopOut = strcmp(NoStopOut,'on');

MaxIter = MaxIterIn*MaxIterOut; % Note: MaxIter can be quite big.

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

if MaxIter>100
    warning(['Due to the choice of K, MaxIterIn and MaxIterOut (by default or by the user),'...
    'IRrestart will potentially perform ', num2str(MaxIter), ' iterations. ',...
    'Please modify K or MaxIterIn and MaxIterOut if you wish to reduce the number of iterations.'])
end

% Checking that the stopping criteria are properly set.
if ischar(NoiseLevel) && strcmp(RegParam,'discrep')
    error('The noise level must be assigned')
end

if ~hybrid && strcmpi(stopOut,'Lxstab')
    warning([' the outer stopping criterion ''Lxstab'' is not meaningful in this case. ',...
            '''Lxstab'' is changed to ''xstab''.'])
    stopOut = 'xstab';
end
if ~hybrid && strcmpi(stopOut,'regPstab')
    warning([' the outer stopping criterion ''regPstab'' can only be used when an inner hybrid method is considered. ',...
            '''regPstab'' is changed to ''xstab''.'])
    stopOut = 'xstab';    
end

% Checking that the constraints are properly set.
if ~( strcmp(adaptConstr, 'nn')||strcmp(adaptConstr, 'box')||strcmp(adaptConstr, 'energy')||strcmp(adaptConstr, 'project') ) && ~strcmp(L, 'identity')
    warning('non-standard regularization terms cannot be handled in this case')
    L = 'identity';
end

if strcmp(adaptConstr, 'box')
    try
        xMin = IRget(options, 'xMin', [], 'fast');
        xMax = IRget(options, 'xMax', [], 'fast');
    catch
        error('When imposing box constraints both a lower bound (''xMin'') and an upper bound (''xMax'') should be provided')
    end
    if strcmp(xMin, 'none') || strcmp(xMax, 'none')
        error('When imposing box constraints both a scalar lower bound (''xMin'') and a scalar upper bound (''xMax'') should be provided')
    end
    if xMin > xMax
        error('The value of ''xMin'' should be smaller or equal than the value of ''xMax''')
    end
elseif strcmp(adaptConstr, 'energy')
    try
        xEnergy = IRget(options, 'xEnergy',   [], 'fast');
    catch
        error('When imposing preservation of volume a value for the volume (''xEnergy'') should be provided')
    end
    if strcmp(xEnergy, 'none')
        error('When imposing preservation of volume a value for the volume (''xEnergy'') should be provided')
    end
elseif strcmp(adaptConstr, 'project')
    try
        xMin    = IRget(options, 'xMin',   [], 'fast');
        xMax    = IRget(options, 'xMax',   [], 'fast');
        xEnergy = IRget(options, 'xEnergy',[], 'fast');
    catch
        error('When projecting ''xMin'', ''xMax'', and ''xEnergy'' should be provided')
    end
    if xMin > xMax
        error('The value of ''xMin'' should be smaller or equal than the value of ''xMax''')
    end
end

if strcmp(inSolver, 'fgmres') && ~strcmp(L, 'identity')
    warning(['IRfgmres cannot handle non-standard regularization terms;\n'...
            'setting the regularization matrix to the identity'])
    L = 'identity';
end

if strcmp(inSolver, 'rrgmres') && ~strcmp(L, 'identity')
    warning(['IRrrgmres cannot handle priorconditioners;\n'...
            'setting the priorconditioning matrix to the identity'])
    L = 'identity';
end

if strcmp(inSolver, 'fgmres')...
        && ( strcmp(adaptConstr, 'tv')||strcmp(adaptConstr, 'tvnn') )
    error('Approximation of the tv regularization term is not possible when the inner solver is fgmres')
end

if (strcmp(inSolver, 'rrgmres')) &&...
        ~((strcmp(adaptConstr, 'nn')||strcmp(adaptConstr, 'box')||strcmp(adaptConstr, 'energy'))||strcmp(adaptConstr, 'project'))
    error('RRGMRES can only handle nonnegativity (''nn''), box constraints (''box''), and volume (''energy'') constraints)')
end

if (strcmp(inSolver, 'cgls')) &&...
        ~((strcmp(adaptConstr, 'sp')||(strcmp(adaptConstr, 'nn')||strcmp(adaptConstr, 'box')||strcmp(adaptConstr, 'energy'))||strcmp(adaptConstr, 'project')))
    error('CGLS can only handle sparsity (''sp''), nonnegativity (''nn''), box constraints (''box''), and volume (''energy'') constraints)')
end

% Checking the dimensions.
if strcmp(inSolver, 'gmres') || strcmp(inSolver, 'fgmres') || strcmp(inSolver, 'rrgmres')
    n = length(b(:));
    test_sq = ones(n,1);
    try
        test_sq = A_times_vec(A, test_sq);
        if (length(test_sq)~=n)
            error('A and b are incopatible; check the size of A and the size of b')
        end
    catch
        error('A and b are incopatible; check the size of A and the size of b')
    end
end

if strcmp(inSolver, 'lsqr') || strcmp(inSolver, 'cgls')
    Atb = Atransp_times_vec(A, b(:)); 
    n = length(Atb(:));
end

if strcmp(adaptConstr, 'spnn') || strcmp(adaptConstr, 'sp')
    if strcmp(SparsTrans, 'none')
        Trans = 1;
    elseif strcmp(SparsTrans, 'dwt')
        wname   = IRget(options, 'wname',   [], 'fast');
        wlevels = IRget(options, 'wlevels', [], 'fast'); 
        if wlevels > 1/2*(log2(n))
            error('The assigned wavelet levels are too high. Make sure that wlevels <= 1/2*(log2(n))')
        end
        Trans = FreqMatrix('dwt', [sqrt(n) sqrt(n)], wname, wlevels);
    end
end

if strcmp(adaptConstr, 'tv') || strcmp(adaptConstr, 'tvnn')
    nsqr = sqrt(n);
    if nsqr ~= floor(nsqr)
        error('Check the size of x: it should be a vectorialization of a square image')
    end
    d = ones(nsqr,1); 
    deriv1 = spdiags([d -d],0:1,nsqr,nsqr);
    dx = kron(deriv1,speye(nsqr));
    dy = kron(speye(nsqr),deriv1);
    D = [dx;dy];
end

% Setting thre regularization matrix (if necessary).
if hybrid 
    if strcmpi(L, 'Laplacian1D')
        L = LaplacianMatrix1D(n);
    elseif strcmpi(L, 'Laplacian2D')
        L = LaplacianMatrix2D(n);
    end
end

x0 = IRget(options, 'x0', [], 'fast');

if strcmp(x0,'none')
    x0 = zeros(n,1);
end
if strcmp(adaptConstr,'nn') || strcmp(adaptConstr,'tvnn') || strcmp(adaptConstr,'spnn')
    x0(x0<0) = 0;
elseif strcmp(adaptConstr,'box')
    x0(x0<xMin) = xMin;
    x0(x0>xMax) = xMax;
elseif strcmp(adaptConstr, 'energy')
    x0(x0<0) = 0;
elseif strcmp(adaptConstr, 'project')
    if xMin<0
        xMin = 0;
        warning('In order to compute the projection, the value of xMin has been changed to 0')
        x0(x0<xMin) = xMin;
    end
end

x0 = x0(:);

% Setting the initial regularization matrix.
if sum(abs(x0)) > 1e-15  
    if strcmp(adaptConstr,'tv') || strcmp(adaptConstr,'tvnn')
        temp_x = dx*x0;
        temp_y = dy*x0;
        diagW  = temp_x.^2+temp_y.^2;
        diagW(diagW<=thr0) = thr0;
        if sum(diagW) <= n*thr0
            L = 'identity';
        else
            Wd=spdiags(1./sqrt(sqrt(diagW)),0:0,n,n);
            W=kron(speye(2),Wd);    
            L=W*D;
        end
    elseif strcmp(adaptConstr,'spnn')
            dp = x0;
            dp(dp < thr0) = thr0;
            dp = 1./sqrt(dp);
            L=spdiags(dp,0:0,n,n);
    elseif strcmp(adaptConstr,'sp')
            dp = abs(x0);
            dp(dp < thr0) = thr0;
            if strcmp(inSolver, 'cgls')
                L = @(x, transp_flag)weightransf(x, dp, Trans, transp_flag);
            else
                dp = 1./sqrt(dp);
                L=spdiags(dp,0:0,n,n);
            end
    elseif strcmp(adaptConstr,'none')
        L = 'identity';
    end
else
    if ~strcmp(adaptConstr, 'nn'), L = 'identity'; end
end

% Declare matrices.
X    = zeros(n,length(K));
Xout = zeros(n,MaxIterOut);
Xnrm = zeros(MaxIter,1);
Rnrm = zeros(MaxIter,1);
if strcmp(inSolver, 'cgls') 
    reqNEres = 1;
    defaultopt.NE_Rtol = 1e-12;
    options = IRset(defaultopt, options);
    NE_Rtol = IRget(options, 'NE_Rtol', [], 'fast');
    NE_Rnrm = zeros(MaxIter,1);
else
    reqNEres = 0;
end
if hybrid, RegP    = zeros(MaxIter,1); end
IterIn  = zeros(MaxIterOut,3);
saved_iterations = zeros(1, length(K));
if strcmp(x_true,'none')
    errornorms = false;
else
    errornorms = true;
    Enrm = zeros(MaxIter,1);
    BestEnrm = 1e10;
end

ktotcount = 0; % Counter for the total iterations
kout = 0; % Counter for the outer iterations
lsi  = 0; % Length of the saved iterations
NotYetStopped    = 1;
StopIt = MaxIter;
StopFlag_out = 'The outer stopping criterion is never satisfied';
StopReg.It = MaxIter;

Kin = min(MaxIterIn, K(end));
% TO BE UPDATED (within cycles)
if  strcmp(warmrestart, 'on')
    paramin.x0 = x0;
else
    paramin.x0 = zeros(n,1);
end
paramin.ktotcount = ktotcount;
% TO BE (most probably) UPDATED (within cycles, depends on the options)
paramin.RegMatrix = L;
paramin.RegParam0 = RegParamk;
% NOT TO BE UPDATED (within cycles)
paramin.Ktot = K;
paramin.RegParam = options.RegParam;
paramin.x_true = options.x_true;
paramin.IterBar = 'off';
paramin.stopGCV = options.stopGCV;
paramin.resflatTol = options.resflatTol;
paramin.GCVflatTol = options.GCVflatTol;
paramin.GCVminTol = options.GCVminTol;
paramin.NoStop = NoStop;
paramin.NoiseLevel = NoiseLevel;
paramin.eta = eta;
paramin.restart = 'on';
paramin.TotIterMax = MaxIter;
paramin.verbosity = verbose;
if (strcmp(inSolver, 'cgls')||strcmp(inSolver, 'rrgmres')) && isscalar(NoiseLevel)
    paramin.Rtol = eta*NoiseLevel;
end
if reqNEres
    paramin.NE_Rtol = NE_Rtol;
end

% Iterate.
noIterBar = strcmp(IterBar, {'off'});
if ~noIterBar
  h_wait = waitbar(0, 'Running iterations, please wait ...');
end

% Potentially iterate till the maximum number of (total) iterations
% is reached.
while ktotcount < MaxIter 
    kout = kout + 1;
    if ~noIterBar
        waitbar(ktotcount/MaxIter, h_wait)
    end
    switch inSolver
        case{'gmres'}
            [Xin, infoin] = IRhybrid_gmres(A, b, Kin, paramin);
        case{'fgmres'}
            [Xin, infoin] = IRhybrid_fgmres(A, b, Kin, paramin);
        case{'lsqr'}
            [Xin, infoin] = IRhybrid_lsqr(A, b, Kin, paramin);
        case{'cgls'}
            [Xin, infoin] = IRcgls(A, b, Kin, paramin);
        case{'rrgmres'}
            [Xin, infoin] = IRrrgmres(A, b, Kin, paramin);
    end
    its_in = infoin.its;
    saved_iterations_in = infoin.saved_iterations;
    lsi_temp = length(saved_iterations_in);
    % imposing constraints on the columns of Xin
    if strcmp(adaptConstr,'nn') || strcmp(adaptConstr,'tvnn') || strcmp(adaptConstr,'spnn')
        Xin(Xin<0) = 0;
    elseif strcmp(adaptConstr,'box')
        Xin(Xin<xMin) = xMin;
        Xin(Xin>xMax) = xMax;
    elseif strcmp(adaptConstr, 'energy')
        Xin(Xin<0) = 0;
        for i = 1:lsi_temp, Xin(:,i) = gdnnf_projection(Xin(:,i), xEnergy); end
    elseif strcmp(adaptConstr, 'project')
        Xin(Xin<xMin) = xMin;
        Xin(Xin>xMax) = xMax;
        for i = 1:lsi_temp, Xin(:,i) = gdnnf_projection(Xin(:,i), xEnergy); end
    end
    % Setting the new x0.
    x0 = Xin(:,lsi_temp);
    Xout(:,kout) = x0;
    x0old = paramin.x0;
    if  strcmp(warmrestart, 'on')
        paramin.x0 = x0;
    else
        paramin.x0 = zeros(n,1);
    end
    test_saved_it = any(saved_iterations_in == K);
    ltest_saved_it = sum(test_saved_it);
    if ltest_saved_it ~= lsi_temp
        saved_iterations_in = saved_iterations_in(test_saved_it);
        lsi_temp = ltest_saved_it;
    end
%     if ~any(saved_iterations_in(lsi_temp) == K)
%         lsi_temp = 0; %lsi_temp-1;
%         saved_iterations_in = saved_iterations_in(1:lsi_temp);
%     end
    X(:,(lsi+1):(lsi+lsi_temp)) = Xin(:,1:lsi_temp);
    saved_iterations((lsi+1):(lsi+lsi_temp)) = saved_iterations_in;
    lsi = lsi + lsi_temp;
    ktot_temp = infoin.ktotcount;
    % Setting the regularization parameter, in the hybrid method case.
    if hybrid
        RegParamkold = RegParamk;
        RegP_in = infoin.RegP;
        RegP(ktotcount+1:ktot_temp) = RegP_in;
        RegParamk = RegP_in(its_in);
        paramin.RegParam0 = RegParamk;
    end
    Xnrm(ktotcount+1:ktot_temp) = infoin.Xnrm;
    Rnrm(ktotcount+1:ktot_temp) = infoin.Rnrm;
    if reqNEres
        NE_Rnrm(ktotcount+1:ktot_temp) = infoin.NE_Rnrm;
    end
    if errornorms
        Enrm(ktotcount+1:ktot_temp) = infoin.Enrm;
        if ~isempty(infoin.BestReg.It)
            BestEnrm_temp = infoin.BestReg.Enrm;
            if BestEnrm_temp < BestEnrm
                BestEnrm = BestEnrm_temp;
                BestReg.It = ktotcount + infoin.BestReg.It;
                if hybrid, BestReg.RegP = infoin.BestReg.RegP; end
                BestReg.X = infoin.BestReg.X;
                BestReg.Enrm = BestEnrm;
            end
        else
            BestReg = [];
        end
    end
    ktotcount = ktot_temp;
    IterIn(kout,:) = [kout, its_in, ktotcount];
    StopFlag_in = infoin.StopFlag;
    if ~strcmp(StopFlag_in,'Reached maximum number of iterations')
        NotYetStopped = 0;
    end
    paramin.ktotcount = ktotcount;
    % Setting the new L.
    Lold = paramin.RegMatrix;
    if strcmp(Lold, 'identity'), Lold = speye(n); end
    singular1 = sum(abs(x0(:)));
    singularTV = 0;
    if singular1 > 1e-15  
        if strcmp(adaptConstr,'tv') || strcmp(adaptConstr,'tvnn')
            temp_x=dx*x0(:);
            temp_y=dy*x0(:);
            diagW=temp_x.^2+temp_y.^2;
            diagW(diagW<=thr0)=thr0;
            if sum(diagW) <= n*thr0
                singularTV = 1;
            end
            Wd=spdiags(1./sqrt(sqrt(diagW)),0:0,n,n);
            W=kron(speye(2),Wd);    
            L=W*D;
            paramin.RegMatrix = L;
        elseif strcmp(adaptConstr,'spnn')
            dp = x0;
            dp(dp < thr0) = thr0;
            dp = 1./sqrt(dp);
            L=spdiags(dp,0:0,n,n);
            paramin.RegMatrix = L;
        elseif strcmp(adaptConstr,'sp')
            dp = abs(x0);
            dp(dp < thr0) = thr0;
            if strcmp(inSolver, 'cgls')
                L = @(x, transp_flag)weightransf(x, dp, Trans, transp_flag);
            else
                dp = 1./sqrt(dp);
                L=spdiags(dp,0:0,n,n);
            end
        end
        paramin.RegMatrix = L;
    end
    if (singular1 <= 1e-15) || singularTV
        disp('The diagonal weighting matrix is numerically zero')
        % X = X(:,1:lsi);
        if lsi == 0
            X = x0;
        else
            X = X(:,1:lsi);
        end
        if nargout==2
            info.its = ktotcount;
            info.Xout = Xout(:,1:kout);
            info.itsInOut = IterIn(1:kout,:);
            % info.saved_iterations = saved_iterations(1:lsi);
            if lsi == 0
                info.saved_iterations = ktotcount;
            else
                info.saved_iterations = saved_iterations(1:lsi);
            end
            info.Rnrm = Rnrm(1:ktotcount);
            info.Xnrm = Xnrm(1:ktotcount);
            if errornorms
                info.Enrm = Enrm(1:ktotcount);
                info.BestReg = BestReg;
            end
            if hybrid, info.RegP = RegP(1:ktotcount); end
            if reqNEres
                info.NE_Rnrm = NE_Rnrm(1:ktotcount);
            end
            if NotYetStopped
                StopFlag_in = 'The stopping criterion of the inner iterations is never satisfied';
            else
                StopFlag_in = 'The stopping criterion is satisfied at least once during the inner iterations';
            end
            info.StopFlag_in = StopFlag_in;
            info.StopFlag_out = 'The diagonal weighting matrix is numerically zero';
            StopReg.It = ktotcount;
            StopReg.X = x0(:); 
            if hybrid, StopReg.RegP = RegParamk; end
            if errornorms
                nrmtrue = norm(x_true(:));
                StopReg.Enrm = norm(x_true(:) - x0(:))/nrmtrue;
            end
            info.StopReg = StopReg;
        end
        if ~noIterBar, close(h_wait), end
        return
    end
    % Check if the stopping criterion is satisfied.
    if StopIt == MaxIter
        if strcmpi(stopOut, 'xstab')
            if norm(x0(:) - x0old(:))/norm(x0old(:)) <= stabOut
                StopIt = ktotcount;
                StopFlag_out = 'The solution stabilizes';
                disp('The solution stabilizes')
            end
        elseif strcmpi(stopOut, 'Lxstab')
            if strcmp(L, 'identity')
                if norm(x0(:) - x0old(:))/norm(x0old(:)) <= stabOut
                    StopIt = ktotcount;
                    StopFlag_out = 'The solution stabilizes';
                    disp('The solution stabilizes')
                end
            else
                Lx0 = L*x0; Lx0old = Lold*x0old;
                if norm(Lx0(:) - Lx0old(:))/norm(Lx0old(:)) < stabOut
                    StopIt = ktotcount;
                    StopFlag_out = 'The transformed solution stabilizes';
                    disp('The transformed solution stabilizes')
                end
            end
        elseif strcmpi(stopOut, 'regPstab')
            if abs(RegParamk - RegParamkold)/RegParamkold < stabOut
                StopIt = ktotcount;
                StopFlag_out = 'The regularization parameter stabilizes';
                disp('The regularization parameter stabilizes')
            end
        end
        if StopIt < MaxIter
            StopReg.It = ktotcount;
            StopReg.X = x0(:); 
            if hybrid, StopReg.RegP = RegParamk; end
            if errornorms
                nrmtrue = norm(x_true(:));
                StopReg.Enrm = norm(x_true(:) - x0(:))/nrmtrue;
            end
            if ~NoStopOut
                % X = X(:,1:lsi);
                if lsi == 0
                    X = x0;
                else
                    X = X(:,1:lsi);
                end
                if nargout==2
                    info.its = ktotcount;
                    info.Xout = Xout(:,1:kout);
                    info.itsInOut = IterIn(1:kout,:);
                    % info.saved_iterations = saved_iterations(1:lsi);
                    if lsi == 0
                        info.saved_iterations = ktotcount;
                    else
                        info.saved_iterations = saved_iterations(1:lsi);
                    end
                    info.Rnrm = Rnrm(1:ktotcount);
                    info.Xnrm = Xnrm(1:ktotcount);
                    if errornorms
                        info.Enrm = Enrm(1:ktotcount);
                        info.BestReg = BestReg;
                    end
                    if hybrid, info.RegP = RegP(1:ktotcount); end
                    if reqNEres
                        info.NE_Rnrm = NE_Rnrm(1:ktotcount);
                    end
                    if NotYetStopped
                        StopFlag_in = 'The stopping criterion of the inner iterations is never satisfied';
                    else
                        StopFlag_in = 'The stopping criterion is satisfied at least once during the inner iterations';
                    end
                    info.StopFlag_in = StopFlag_in;
                    info.StopFlag_out = StopFlag_out;
                    info.StopReg = StopReg;
                end
                if ~noIterBar, close(h_wait), end
                return
            end
        end
    end  
end
% Finalising the iterations.
IterIn = IterIn(1:kout, :);
if ~noIterBar
    waitbar(ktotcount/MaxIter, h_wait);
    close(h_wait)
end
% X = X(:, 1:lsi);
if lsi == 0
    X = x0;
else
    X = X(:,1:lsi);
end
% saved_iterations = saved_iterations(1:lsi);
if lsi == 0
    info.saved_iterations = ktotcount;
else
    info.saved_iterations = saved_iterations(1:lsi);
end
Xout = Xout(:, 1:kout);
if StopIt == MaxIter
    StopReg.It = ktotcount;
    StopReg.X = x0(:); 
    if hybrid, StopReg.RegP = RegParamk; end
    if errornorms
        nrmtrue = norm(x_true(:));
        StopReg.Enrm = norm(x_true(:) - x0(:))/nrmtrue;
    end
end    
if nargout==2
  info.its = ktotcount;
  info.saved_iterations = saved_iterations;
  info.itsInOut = IterIn;
  info.Xout = Xout;
  if NotYetStopped
    StopFlag_in = 'The stopping criterion of the inner iterations is never satisfied';
  else
    StopFlag_in = 'The stopping criterion is satisfied at least once during the inner iterations';
  end
  info.StopFlag_in = StopFlag_in;
  info.StopFlag_out = StopFlag_out;
  info.StopReg = StopReg;
  info.Rnrm = Rnrm(1:ktotcount);
  info.Xnrm = Xnrm(1:ktotcount);
  if errornorms
    info.Enrm = Enrm(1:ktotcount);
    info.BestReg = BestReg;
  end
  if hybrid, info.RegP = RegP(1:ktotcount); end
  if reqNEres
      info.NE_Rnrm = NE_Rnrm(1:ktotcount);
  end
end