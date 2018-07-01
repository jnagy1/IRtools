function [X, info] = IRhybrid_fgmres(A, b, varargin)
%IRhybrid_fgmres Hybrid version of FGMRES for enforcing 1-norm penalization
%
% options  = IRhybrid_fgmres('defaults')
% [X,info] = IRhybrid_fgmres(A,b)
% [X,info] = IRhybrid_fgmres(A,b,K)
% [X,info] = IRhybrid_fgmres(A,b,options)
% [X,info] = IRhybrid_fgmres(A,b,K,options)
%
% IRhybrid_fgmres is a hybrid iterative regularization method used for 
% enforcing sparsity on the solution of large-scale, ill-posed inverse 
% problems of the form:
%               b = A*x + noise,    A square matrix.
% Sparsity is achieved by including a 1-norm penality term, which is 
% implemented via an iteratively reweigthed 2-norm penality term. Due to
% a change of variable, the weighting matrix is regarded as an iteration-
% dependent right preconditioner for A. The method combines a flexible
% GMRES iteration with a Tikhonov regularization method.
%
% With 'defaults' as input returns the default options.  Otherwise outputs
% the iterates specified in K, using max(K) as MaxIter, and using all other
% default options.  With options as input: uses the user-specified options
% and all the other default options.
%
% Inputs:
%  A : either (a) a full or sparse matrix (must be square)
%             (b) a matrix object that performs the matrix*vector operation
%             (c) user-defined function handle
%  b : right-hand side vector
%  K : (optional) integer vector that specifies which iterates are returned
%      in X; the maximum number of iterations is assumed to be max(K)
%      [ positive integer | vector of positive components ]
%  options : structure with the following fields (optional)
%      x0         - initial guess for the iterations; default = zero vector
%                   [ array | {'none'} ]
%      MaxIter    - maximum allowed number of cgls iterations
%                   [ positive integer | {100} ]
%                   NOTE: K overrules MaxIter if both are assigned
%      x_true     - true solution; allows us to returns error norms
%                   with respect to x_true at each iteration
%                   [ array | {'none'} ]
%      RegParam   - a value or a method to find the regularization
%                   parameter for the projected problems: 
%                   [non-negative scalar | {'gcv'} | 'discrep' ]
%                   This also determines which stopping rule is used
%                   If 'gcv' is chosen, the iteration is stopped when the
%                     GCV function minimum stabilizes or increases within
%                     a certain window of iterations (see 'stopGCV',
%                     'FlatTol' and 'MinTol').
%                   If 'discrep' is chosen, and NoiseLevel is provided then
%                     the discrepancy principle is as stopping criterion
%                     (see 'NoiseLevel' and 'eta').
%      stopGCV    - stopping criterion for the iterations when GCV is used
%                   [ 'GCVvalues' | {'resflat'} ]
%      FlatTol    - tolerance for detecting flatness (stabilization)
%                   in the GCV function as a stopping criterion
%                   [ {10^-6} | non-negative scalar ]
%      MinTol     - window of iterations: if the GCV minimum continues to
%                   increase over this window, then the iterations are
%                   stopped:
%                   [ {3} | positive integer ]
%      NoiseLevel - norm of noise in rhs divided by norm of rhs (must be
%                   assigned if RegParam is 'discrep')
%                   [ {'none'} | nonnegative scalar ]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      RegParam0  - regularization parameter used in the first 
%                   projected problem (needed if RegParam is 'discrep')
%                   [ {1} | positive scalar ]
%      tolX       - tolerance for the numerically zero components in the
%                   preconditioning matrix
%                   [ {10^-10} | non-negative scalar ]
%      DecompOut  - returns the flexible Arnoldi decomposition to the user
%                   [ 'on' | {'off'} ]
%      IterBar    - shows the progress of the iterations
%                   [ {'on'} | 'off' ]
%      NoStop     - specifies whether the iterations should proceed
%                   after a stopping criterion is satisfied
%                   [ 'on' | {'off'} ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its      - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - string that describes the output/stopping condition:
%                   * Flat GCV curve 
%                   * Min of GCV function (within window of MinTol its)
%                   * Performed max number of iterations
%                   * Discrepancy principle satisfied
%                   * Breakdown of the Arnoldi algorithm
%                   * The residual norm stabilizes
%                   * Reached maximum number of iterations
%      Rnrm     - relative residual norms at each iteration
%      Xnrm     - solution norms
%      Enrm     - relative error norms (requires x_true) at each iteration
%      StopIt   - number of the iteration satisfying the stopping criterion
%      RegP     - sequence of the regularization parameters
%      Z        - Flexible Arnoldi basis vectors
%      V        - Orthonomal vectors in range(A)
%      H        - upper Hessenberg matrix from Arnoldi
%      GCValues - GCV function evaluated at the minimum at each iteration
%
% See also: IRell1, IRhybrid_gmres, IRhybrid_lsqr, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('RegParam', 'gcv', 'x0', 'none', 'x_true', 'none', ...
    'MaxIter', 100 , 'IterBar', 'on', 'stopGCV', 'resflat', ...
    'resflatTol', 0.05, 'GCVflatTol', 10^-6, 'GCVminTol', 3, ...
    'tolX', 10^-10, 'NoStop', 'off', 'NoiseLevel', 'none', 'eta', 1.01, ...
    'RegParam0', 1, 'DecompOut', 'off');

if nargin == 0
    error('Not enough input arguments')
elseif nargin == 1 
    % If input is 'defaults,' return the default options in X.
    if nargout <= 1 && isequal(A,'defaults')
        X = defaultopt;
        return;
    else
        error('Not enough input arguments')
    end
end

% This function includes an undocumented restart feature which allows the
% functions IRrestart to use IRhybrid_fgmres in a restarted fashion.
defaultopt.restart = 'off';
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

RegParam  = IRget(options, 'RegParam',  [], 'fast');
MaxIter   = IRget(options, 'MaxIter',   [], 'fast');
x_true    = IRget(options, 'x_true',    [], 'fast');
IterBar   = IRget(options, 'IterBar',   [], 'fast');
stopGCV   = IRget(options, 'stopGCV',   [], 'fast');
NoiseLevel= IRget(options, 'NoiseLevel',[], 'fast');
eta       = IRget(options, 'eta',       [], 'fast');
RegParamk = IRget(options, 'RegParam0', [], 'fast');
resdegflat= IRget(options, 'resflatTol',[], 'fast');
degflat   = IRget(options, 'GCVflatTol',[], 'fast');
tolX      = IRget(options, 'tolX',      [], 'fast');
mintol    = IRget(options, 'GCVminTol', [], 'fast');
NoStop    = IRget(options, 'NoStop',    [], 'fast');
restart   = IRget(options, 'restart',   [], 'fast');
verbose   = IRget(options, 'verbosity', [], 'fast');
DecompOut = IRget(options, 'DecompOut', [], 'fast');

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

if ischar(NoiseLevel) && strcmp(RegParam,'discrep')
    error('The noise level must be assigned')
end

StopIt = MaxIter;

restart = strcmp(restart, 'on');
if restart
    ktotcount  = IRget(options, 'ktotcount', [],'fast');
    TotIterMax = IRget(options, 'TotIterMax',[],'fast');
    if strcmp(TotIterMax, 'none') || TotIterMax < MaxIter
        TotIterMax = MaxIter;
    end
    if strcmp(ktotcount, 'none')
        error('the total iteration counter must be assigned')
    end
    Ktot = IRget(options, 'Ktot', [],'fast');
end

n = length(b(:));
test_sq = ones(n,1);
try
    test_sq = A_times_vec(A, test_sq);
    if (length(test_sq)~=n)
        error('The matrix A shuold be square; check the length of b.')
    end
catch
    error('The matrix A must be square; check the length of b.')
end

% Setting x0.
x0 = IRget(options, 'x0', [], 'fast');
if strcmp(x0,'none')
    r = b(:);
    x0 = zeros(n,1);
    precX = ones(n,1);
else
    if length(x0(:))~=n
        error('Check the length of x0')
    end
    if max(abs(x0(:))) == 0
        r = b;
        x0 = zeros(n,1);
        precX = ones(n,1);
    else
        try
            Ax0 = A_times_vec(A, x0);
            if (length(test_sq)~=n)
                error('The matrix A shuold be square; check the length of b.')
            end
            r = b - Ax0;
            precX = x0;
            precX(precX < tolX) = eps;
            precX = sqrt(precX);
        catch
            error('The matrix A must be square; check the length of b.')
        end    
    end
end
beta = norm(r);
nrmb = norm(b(:));

% There is no true solution.
notrue = strcmp(x_true,'none');
% We do not want to stop when the stopping criterion is satisfied.
NoStop = strcmp(NoStop,'on');

% Assessing if we want inner Tikhonov regularization.
if strcmp(RegParam,'off')
    RegParam = 0;
end

% Declare matrices.
X                = zeros(n,length(K));
Xnrm             = zeros(max(K),1);
Rnrm             = zeros(max(K),1);
RegParamVect     = zeros(max(K),1);
h                = zeros(max(K),1); % New column of Hessenberg matrix H
H                = zeros(max(K)+1,max(K)); 
Z                = zeros(n,max(K));
V                = zeros(n,max(K)+1);
rhs              = zeros(max(K)+1,1); % Projected right-hand side
if restart
    saved_iterations = zeros(1, length(Ktot));
else
    saved_iterations = zeros(1, length(K));
end
GCV = zeros(max(K),1);
warningGCV = 0;
if notrue
    errornorms = false;
else
    errornorms = true;
    nrmtrue = norm(x_true(:));
    Enrm = zeros(max(K),1);
    BestReg.RegP = [];
    BestReg.It = [];
    BestReg.X =[];
    BestReg.Enrm = [];
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestEnrm = 1e10;
end

% Main code begins here.
V(:,1) = r/beta;
rhs(1) = beta;

% Iterate.
noIterBar = strcmp(IterBar,{'off'});
if ~noIterBar
  h_wait = waitbar(0, 'Running iterations, please wait ...');
end
j = 0;
for k=1:MaxIter
    if ~noIterBar
        waitbar(k/MaxIter, h_wait)
    end 
    if restart, ktotcount = ktotcount + 1; end
    % Flexible Arnoldi algorithm.
    Z(:,k) = precX.*V(:,k);
    w = Z(:,k);    
    w = A_times_vec(A, w);
    % Modified Gram-Schmidt on the new vector.
    for ll=1:k
        h(ll) = V(:,ll)'*w;
        w = w - V(:,ll)*h(ll);
    end
    alpha = norm(w);
    % Store new Arnoldi vector and update projected rhs.
    V(:,k+1) = w/alpha;
    rhsk = rhs(1:k+1); % Current projected rhs.
    H(1:k+1,k) = [h(1:k); alpha];    
    if abs(alpha) <= eps
        if verbose
            disp('Hessenberg matrix is (numerically) singular')
        end
        H = H(1:k+1,1:k);
        Z = Z(:,1:k);
        V = V(:,1:k);
        X(:,j+1) = x;
        X = X(:,1:j+1);
        if restart
            saved_iterations(j+1) = ktotcount-1;
        else
            saved_iterations(j+1) = k-1;
        end
        saved_iterations = saved_iterations(1:j+1);
        if k-1 < StopIt, StopIt = k-1; end
        if k>1
            Xnrm    = Xnrm(1:k-1);
            Rnrm    = Rnrm(1:k-1);
            RegParamVect    = RegParamVect(1:k-1);
            if errornorms, Enrm = Enrm(1:k-1); end
        end
        % Stop because the Hessenberg matrix is (numerically) rank deficient.
        % No chioce: even if NoStop is 'on'...we simpy cannot compute the
        % solution anymore.
        if StopIt == MaxIter
            StopFlag = 'Breakdown of the Arnoldi algorithm';
            StopReg.RegP = RegParamk;
            StopReg.It = k-1;
            StopReg.X = x;
            if errornorms, StopReg.Enrm = Enrm(k-1); end
        end
        break
    end
    Hk = H(1:k+1,1:k);
    [Uk, Sk, Vk] = svd(Hk);
    if k==1
        Sk = Sk(1,1);
    else
        Sk = diag(Sk);
    end
    rhskhat = Uk'*rhsk;
    gmres_res = abs(Uk(:,k+1)'*rhsk)/nrmb;
    if isscalar(RegParam)
        RegParamk = RegParam;
        RegParamVect(k) = RegParamk;
    elseif strcmp(RegParam,'discrep')
        if k==1 
            RegParamVect(k) = RegParamk;
        end
    elseif strcmp(RegParam,'gcv')
        RegParamk = fminbnd('TikGCV', 0, Sk(1), [], rhskhat, Sk, n);
        RegParamVect(k) = RegParamk;
    else
        error('Invalid parameter choice method')
    end
    Dk = Sk.^2 + RegParamk^2;
    rhskhat = Sk .* rhskhat(1:k);
    yhat = rhskhat ./ Dk;
    y = Vk * yhat;
    Rnrm(k) = norm(rhsk - Hk*y)/nrmb;
    d = Z(:,1:k)*y;
    x = x0 + d;
    % Compute norms.
    Xnrm(k) = norm(x(:));
    if errornorms
        Enrm(k) = norm(x_true-x)/nrmtrue;
        if Enrm(k)<BestEnrm
            BestReg.RegP = RegParamk;
            BestReg.It = k;
            BestReg.X = x;
            BestEnrm = Enrm(k);
            BestReg.Enrm = BestEnrm;
            BestReg.Xnrm = Xnrm(k);
            BestReg.Rnrm = Rnrm(k);
        end
    end
    precX = abs(x);
    precX(precX < tolX) = eps;
    precX = sqrt(precX);
    AlreadySaved = 0;
    if any(k==K)
        j = j+1;
        X(:,j) = x;  
        saved_iterations(j) = k;
        % This is used to save the last iteration, in case K = MaxIterIn
        % (when performing restarts, and the inner stopping criterion is
        % not satisfied)
        if restart, saved_iterations(j) = ktotcount; end
        AlreadySaved = 1;              
    end
    if restart
        if any(ktotcount == Ktot) && ~ AlreadySaved
            saved_iterations(j) = ktotcount;
            AlreadySaved = 1;                
        end
        if ktotcount == TotIterMax
            if ~ AlreadySaved
                j = j+1;
                saved_iterations(j) = ktotcount;
                X(:,j) = x; 
            end
            StopIt = k;
            StopReg.RegP = RegParamk;
            StopReg.It = k;
            StopReg.X = x;  
            if errornorms
                StopReg.Enrm = Enrm(k);
                Enrm = Enrm(1:k);
            end
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            RegParamVect    = RegParamVect(1:k);
            H = H(1:k+1,1:k);
            X = X(:,1:j);
            saved_iterations = saved_iterations(1:j);
            if verbose
                disp('reached maximum number of iterations')
            end
            StopFlag = 'reached maximum number of iterations';
            break
        end
    end 
    if isscalar(RegParam)
        % Purely iterative method case.
        if strcmp(NoiseLevel, 'none')
            if k>1
            if abs((Rnrm(k)-Rnrm(k-1)))/Rnrm(k-1) < resdegflat && ...
                Rnrm(k) == min(Rnrm(1:k)) && StopIt == MaxIter
                if verbose
                    disp('The stopping criterion for fgmres is satisfied')
                end
                % Stop because the residual stabilizes.
                StopFlag = 'The residual norm stabilizes';
                if ~AlreadySaved && ~NoStop
                    j = j+1;
                    X(:,j) = x;
                    if restart
                        saved_iterations(j) = ktotcount;
                    else
                        saved_iterations(j) = k;
                    end
                    AlreadySaved = 1;
                end
                StopIt = k;
                StopReg.RegP = RegParamk;
                StopReg.It = k;
                StopReg.X = x;
                if errornorms, StopReg.Enrm = Enrm(k); end
                if ~ NoStop
                    Xnrm    = Xnrm(1:k);
                    Rnrm    = Rnrm(1:k);
                    RegParamVect = RegParamVect(1:k);
                    H = H(1:k+1,1:k);
                    Z = Z(:,1:k);
                    V = V(:,1:k+1);
                    if errornorms, Enrm = Enrm(1:k); end
                    X = X(:,1:j);
                    saved_iterations = saved_iterations(1:j);
                    break
                end
            end
            end
        else
            if Rnrm(k) < eta*NoiseLevel
            % Stopping criterion.
            if StopIt == MaxIter
                if verbose
                    disp('The discrepancy principle is satisfied')
                end
                StopFlag = 'The discrepancy principle satisfied';
                if ~AlreadySaved && ~NoStop
                    j = j+1;
                    X(:,j) = x;
                    if restart
                        saved_iterations(j) = ktotcount;
                    else
                        saved_iterations(j) = k;
                    end
                    AlreadySaved = 1;
                end
                StopIt = k;
                StopReg.RegP = RegParamk;
                StopReg.It = k;
                StopReg.X = x;
                if errornorms, StopReg.Enrm = Enrm(k); end
                if ~ NoStop
                    Xnrm    = Xnrm(1:k);
                    Rnrm    = Rnrm(1:k);
                    RegParamVect    = RegParamVect(1:k);
                    H = H(1:k+1,1:k);
                    Z = Z(:,1:k);
                    V = V(:,1:k+1);
                    if errornorms, Enrm = Enrm(1:k); end
                    X = X(:,1:j);
                    saved_iterations = saved_iterations(1:j);
                    % Stop because the discrepancy principle is satisfied.
                    break
                end
            end
            end
        end
    elseif strcmp(RegParam,'discrep')
        if Rnrm(k) < eta*NoiseLevel
                % Stopping criterion.
                if StopIt == MaxIter % The method has not stopped yet.
                    if verbose
                        disp('The discrepancy principle is satisfied')
                    end
                    StopFlag = 'The discrepancy principle satisfied';
                    if ~AlreadySaved && ~NoStop
                        j = j+1;
                        X(:,j) = x;
                        if restart
                            saved_iterations(j) = ktotcount;
                        else
                            saved_iterations(j) = k;
                        end
                        AlreadySaved = 1;
                    end
                    StopIt = k;
                    StopReg.RegP = RegParamk;
                    StopReg.It = k;
                    StopReg.X = x;
                    if errornorms, StopReg.Enrm = Enrm(k); end
                    if ~ NoStop
                        Xnrm    = Xnrm(1:k);
                        Rnrm    = Rnrm(1:k);
                        RegParamVect    = RegParamVect(1:k);
                        H = H(1:k+1,1:k);
                        Z = Z(:,1:k);
                        V = V(:,1:k+1);
                        if errornorms, Enrm = Enrm(1:k); end
                        X = X(:,1:j);
                        saved_iterations = saved_iterations(1:j);
                        % Stop because the discrepancy principle is satisfied.
                        break
                    else
                        RegParamk = abs((eta*NoiseLevel - gmres_res)/(Rnrm(k) - gmres_res))*(RegParamk^2);
                        RegParamk = sqrt(RegParamk);
                        if k~=MaxIter, RegParamVect(k+1) = RegParamk; end
                    end
                else
                    RegParamk = abs((eta*NoiseLevel - gmres_res)/(Rnrm(k) - gmres_res))*(RegParamk^2);
                    RegParamk = sqrt(RegParamk);
                    if k~=MaxIter, RegParamVect(k+1) = RegParamk; end
                end
            else
                RegParamk = abs((eta*NoiseLevel - gmres_res)/(Rnrm(k) - gmres_res))*(RegParamk^2);
                RegParamk = sqrt(RegParamk);
                if k~=MaxIter, RegParamVect(k+1) = RegParamk; end
        end
    elseif strcmp(RegParam,'gcv')
        % Check the stopping criterion (all the possibilities).
        GCV(k) = GCVstopfun(RegParamk, Uk(1,:)', Sk, beta, n, n);
        if k > 1
            if strcmpi(stopGCV, 'GCVvalues')
            if abs((GCV(k)-GCV(k-1)))/GCV(1) < degflat && StopIt == MaxIter
                if verbose
                    disp('The stopping criterion for GCV principle is satisfied')
                end
                % Stop because the GCV function is too flat.
                StopFlag = 'GCV function too flat';
                if ~AlreadySaved && ~NoStop
                    j = j+1;
                    X(:,j) = x;
                    if restart
                        saved_iterations(j) = ktotcount;
                    else
                        saved_iterations(j) = k;
                    end
                    AlreadySaved = 1;
                end
                StopIt = k;
                StopReg.RegP = RegParamk;
                StopReg.It = k;
                StopReg.X = x;                
                if errornorms, StopReg.Enrm = Enrm(k); end
                if ~ NoStop
                    Xnrm    = Xnrm(1:k);
                    Rnrm    = Rnrm(1:k);
                    RegParamVect    = RegParamVect(1:k);
                    H = H(1:k+1,1:k);
                    Z = Z(:,1:k);
                    V = V(:,1:k+1);
                    if errornorms, Enrm = Enrm(1:k); end
                    X = X(:,1:j);
                    saved_iterations = saved_iterations(1:j);
                    % Stop because the GCV stopping criterion is satisfied.
                    break
                end
            elseif GCV(k-1) < GCV(k) && ~ warningGCV && StopIt == MaxIter
                % Potential minimum reached. 
                warningGCV = 1;
                % Save data just in case.
                x_save = x;
                k_save = k; % For computing the GCV stopping criterion.
                RegParamk_save = RegParamk;
                if restart, ktotcount_save = ktotcount; end
            elseif warningGCV && length(GCV) > k_save + mintol && StopIt == MaxIter % Passed window
                if GCV(k_save) < GCV(k_save+1:end)
                if verbose
                    disp('The stopping criterion for GCV principle is satisfied')
                end
                StopFlag = 'Increasing GCV minima';
                if ~AlreadySaved && ~NoStop
                    j = j+1;
                    X(:,j) = x_save;
                    if restart
                        saved_iterations(j) = ktotcount_save;
                    else
                        saved_iterations(j) = k_save;
                    end
                    AlreadySaved = 1;
                end
                StopIt = k_save;
                StopReg.RegP = RegParamk_save;
                StopReg.It = k_save;
                StopReg.X = x_save;
                if errornorms, StopReg.Enrm = Enrm(k_save); end
                if ~ NoStop
                    Xnrm    = Xnrm(1:k_save);
                    Rnrm    = Rnrm(1:k_save);
                    RegParamVect    = RegParamVect(1:k_save);
                    H = H(1:k_save+1,1:k_save);
                    Z = Z(:,1:k_save);
                    V = V(:,1:k_save+1);
                    if errornorms, Enrm = Enrm(1:k_save); end
                    X = X(:,1:j);
                    saved_iterations = saved_iterations(1:j);
                    % Stop because the GCV stopping criterion is satisfied.
                    break
                end     
                end
            end
            elseif strcmpi(stopGCV, 'resflat')
            if abs((Rnrm(k)-Rnrm(k-1)))/Rnrm(k-1) < resdegflat && ...
                Rnrm(k) == min(Rnrm(1:k)) && StopIt == MaxIter
                if verbose
                    disp('The stopping criterion for GCV principle is satisfied')
                end
                % Stop because the discrepancy (i.e., the residual for the
                % regularized problem) stabilizes.
                StopFlag = 'the discrepancy norm stabilizes';
                if ~AlreadySaved && ~NoStop
                    j = j+1;
                    X(:,j) = x;
                    if restart 
                        saved_iterations(j) = ktotcount;
                    else
                        saved_iterations(j) = k;
                    end
                    AlreadySaved = 1;
                end
                StopIt = k;
                StopReg.RegP = RegParamk;
                StopReg.It = k;
                StopReg.X = x;
                if errornorms, StopReg.Enrm = Enrm(k); end
                if ~ NoStop
                    Xnrm    = Xnrm(1:k);
                    Rnrm    = Rnrm(1:k);
                    RegParamVect    = RegParamVect(1:k);
                    H = H(1:k+1,1:k);
                    Z = Z(:,1:k);
                    V = V(:,1:k+1);
                    if errornorms, Enrm = Enrm(1:k); end
                    X = X(:,1:j);
                    saved_iterations = saved_iterations(1:j);
                    break
                end
            end
            end
        end            
    end
end

if k == MaxIter 
    if StopIt == MaxIter
        % Stop because max number of iterations reached.
        if verbose
        disp('Reached maximum number of iterations')
        end
        StopFlag = 'reached maximum number of iterations';
        if ~AlreadySaved
            j = j+1;
            X(:,j) = x;
            if restart
                saved_iterations(j) = ktotcount;
            else
                saved_iterations(j) = k;
            end
        end
        StopReg.RegP = RegParamk;
        StopReg.It = k;
        StopReg.X = x;
        if errornorms, StopReg.Enrm = Enrm(k); end
        Xnrm    = Xnrm(1:k);
        Rnrm    = Rnrm(1:k);
        RegParamVect    = RegParamVect(1:k);
        H = H(1:k+1,1:k);
        Z = Z(:,1:k);
        V = V(:,1:k+1);
        if errornorms, Enrm = Enrm(1:k); end
        X = X(:,1:j);
        saved_iterations = saved_iterations(1:j);
    end 
end

if ~noIterBar, close(h_wait), end

if nargout==2
  info.its = k;
  info.StopReg = StopReg;
  info.StopFlag = StopFlag;
  info.Rnrm = Rnrm;
  info.Xnrm = Xnrm;
  info.RegP = RegParamVect;
  if errornorms
    info.Enrm = Enrm;
    info.BestReg = BestReg;
  end
  info.saved_iterations = saved_iterations(1:j);
  if strcmp(RegParam,'gcv')
    info.GCValues = GCV;
  end
  if strcmp(DecompOut,'on')
      info.Z = Z;
      info.V = V;
      info.H = H;
  end
  if restart
      info.ktotcount = ktotcount;
  end
end
