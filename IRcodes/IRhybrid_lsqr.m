function [X, info] = IRhybrid_lsqr(A, b, varargin)
%IRhybrid_lsqr Hybrid version of LSQR algorithm
%
% options  = IRhybrid_lsqr('defaults')
% [X,info] = IRhybrid_lsqr(A,b)
% [X,info] = IRhybrid_lsqr(A,b,K)
% [X,info] = IRhybrid_lsqr(A,b,options)
% [X,info] = IRhybrid_lsqr(A,b,K,options)
%
% IRhybrid_lsqr is a hybrid iterative regularization method used for 
% solving large-scale, ill-posed inverse problems of the form:
%               b = A*x + noise .
% The method combines LSQR iteration (iterative regularization method) 
% with a Tikhonov regularization method to stabilize the semiconvergence
% behavior that is characteristic of many iterative solvers applied to
% ill-posed problems.
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
%      x_true     - true solution; allows us to returns error norms with
%                   respect to x_true at each iteration
%                   [ array | {'none'} ]
%      RegParam   - a value or a method to find the regularization
%                   parameter for the projected problems: 
%                   [  non-neg. scalar | {'wgcv'} | 'gcv' | 'modgcv' | 'discrep' ]
%                   This also determines which stopping rule is used
%                   If 'gcv', 'wgcv' or 'modgcv' is chosen, the iteration is
%                     stopped when the GCV function minimum stabilizes or
%                     increases within a certain window of iterations (see
%                     'stopGCV', 'GCVflatTol' and 'GCVminTol').
%                   If 'discrep' is chosen, and NoiseLevel is provided,
%                     then the discrepancy principle is used as stopping
%                     (see 'NoiseLevel' and 'eta').
%      stopGCV    - stopping criterion for the iterations when GCV is used
%                   [ {'GCVvalues'} | 'resflat' ]
%      resflatTol - tolerace for the stabilization of the residual
%                   (to be used if stopGCV is 'resflat')
%                   [ {0.05} | non-negative scalar ]
%      GCVflatTol - tolerance for detecting flatness (stabilization)
%                   in the GCV function as a stopping criterion
%                   [ {10^-6} | non-negative scalar ]
%      GCVminTol  - window of iterations - if the GCV minimum continues
%                   to increase over this window, then the iterations are
%                   stopped:
%                   [ {3} | positive integer ]
%      GCVweight  - weight to be used if RegParam is 'wgcv'
%                   [ 'adapt' | non-negative scalar ]
%      NoiseLevel - norm of noise in rhs divided by norm of rhs (must be
%                   assigned if RegParam is 'discrep')
%                   [ {'none'} | nonnegative scalar ]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      RegParam0  - regularization parameter used in the first  projected
%                   problem (needed if RegParam is 'discrep')
%                   [ {1} | positive scalar ]
%      MaxIter    - maximum number of iterations
%                   [ {'none'} | positive integer ]
%      RegMatrix  - priorconditioner for the iterations
%                   [ {'identity'} | 'Laplacian1D' | 'Laplacian2D' |
%                   square nonsingular matrix | function handle ]
%      DecompOut  - returns the Golub-Kahan decomposition to the user
%                   [ 'on' | {'off'} ]
%      Reorth     - indicates if reorthogonalization should be
%                   applied to the Golub-Kahan bidiagonalization algorithm
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
%                   * Minimum of GCV function (within window of MinTol its)
%                   * Performed max number of iterations
%                   * Discrepancy principle satisfied
%                   * Breakdown of the Golub-Kahan bidiagonalization algorithm
%      StopReg  - structure with the following fields:
%                   * X: solution satisfying the stopping criterion
%                   * It: iteration satisfying the stopping criterion
%                   * RegP: regularization parameter at the iteration satisfying 
%                     the stopping crierion
%                   * Xnrm: norm of the solution satisfying satisfying the
%                     stopping criterion 
%                   * Rnrm: relative residual norm at the iteration
%                     satisfying the stopping criterion
%                   * Enrm: relative error norm at the iteration
%                     satisfying the stopping criterion (requires x_true)
%      Xnrm     - solution norms at each iteration
%      Rnrm     - relative residual norms at each iteration
%      Enrm     - relative error norms at each iteration (requires x_true)
%      RegP     - sequence of the regularization parameters
%      GCValues - GCV function evaluated at minimum point at each
%                 iteration (if RegParam is 'gcv', 'wgcv', 'modgcv')
%      V        - Golub-Kahan bidiagonalization basis vectors for the solution
%      U        - Golub-Kahan bidiagonalization basis vectors
%      B        - lower bidiagonal matrix computed by Golub-Kahan bidiagonalization
%
% See also: IRcgls, IRhybrid_fgmres, IRhybrid_gmres, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Initialization
defaultopt = struct('x0', 'none', 'MaxIter', 100 ,...
    'x_true', 'none', 'NoStop','off', 'IterBar', 'on',...
    'Reorth', 'off', 'RegParam','wgcv',...
    'RegMatrix', 'Identity', 'GCVweight', 'adapt',...
    'GCVflatTol', 10^-6, 'GCVminTol', 3,...
    'stopGCV', 'GCVvalues', 'resflatTol', 0.05,...
    'NoiseLevel', 'none', 'eta', 1.01, 'RegParam0', 1, 'DecompOut', 'off');

if nargin == 0
    error('Not enough input arguments')
elseif nargin == 1 
    % If input is 'defaults,' return the default options in X
    if nargout <= 1 && isequal(A,'defaults')
        X = defaultopt;
        return;
    else
        error('Not enough input arguments')
    end
end

defaultopt.restart = 'off';
defaultopt.verbosity = 'on';

% Check for acceptable number of optional input arguments
switch length(varargin)
    case 0 
        K = []; options = [];
    case 1
        if isa(varargin{1}, 'double')
            K = varargin{1}; options = [];
        else
            % no matter the order of appearance
            K = []; options = varargin{1};
        end
    case 2
        if isa(varargin{1}, 'double')
            K = varargin{1}; options = varargin{2};
        else
            % again, no matter the order of appearance
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
RegParam   = IRget(options, 'RegParam',   [], 'fast');
x_true     = IRget(options, 'x_true',     [], 'fast');
NoStop     = IRget(options, 'NoStop',     [], 'fast');
Reorth     = IRget(options, 'Reorth',     [], 'fast');
IterBar    = IRget(options, 'IterBar',    [], 'fast');
L          = IRget(options, 'RegMatrix',  [], 'fast');
omega      = IRget(options, 'GCVweight',  [], 'fast');
stopGCV    = IRget(options, 'stopGCV',    [], 'fast');
resdegflat = IRget(options, 'resflatTol', [], 'fast');
degflat    = IRget(options, 'GCVflatTol', [], 'fast');
mintol     = IRget(options, 'GCVminTol',  [], 'fast');
NoiseLevel = IRget(options, 'NoiseLevel', [], 'fast');
eta        = IRget(options, 'eta',        [], 'fast');
RegParamk  = IRget(options, 'RegParam0',  [], 'fast');
restart    = IRget(options, 'restart',    [], 'fast');
verbose    = IRget(options, 'verbosity',  [], 'fast');
DecompOut  = IRget(options, 'DecompOut', [], 'fast');

verbose = strcmp(verbose, 'on');

adaptWGCV = strcmp(RegParam, {'wgcv'}) && strcmp(omega, {'adapt'});

% setting K
if isempty(K)
    K = MaxIter;
end
% sorting the iterations (in case they are shuffled in input)
K = K(:); K = sort(K,'ascend'); K = unique(K);
if ~((isreal(K) && (all(K > 0)) && all(K == floor(K))))
    error('K must be a vector of positive real integers')
end
if K(end) ~= MaxIter
    MaxIter = K(end);    
end
% note that there is no control on K, as it does not go through IRset

if strcmp(Reorth,'on')
    reorth = true; 
else
    reorth = false;
end

if strcmp(RegParam,'discrep') && ischar(NoiseLevel)
    error('The noise level (NoiseLevel) must be assigned')
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
    Ktot = IRget(options, 'Ktot', [], 'fast');
    % no checks on Ktot, it should be given from IRrestart
end

d = Atransp_times_vec(A, b);
n = length(d);
m = length(b);

% setting x0
x0 = IRget(options, 'x0', [], 'fast');
if strcmp(x0,'none')
    x0 = zeros(n,1);
    r = b(:); %% check!!
else
    try
        Ax0 = A_times_vec(A, x0);
        r = b(:) - Ax0;
    catch
        error('Check the length of x0')
    end
end
x = x0;
beta = norm(r(:)); %% check!!
nrmb = norm(b(:)); %% check!!

% means no true solution
notrue = strcmp(x_true,'none');
% means we do not want to stop when the stopping criterion is satisfied
NoStop = strcmp(NoStop,'on');

% assessing if we want inner Tikhonov regularization
if strcmp(RegParam,'off')
    tik = false;
else
    tik = true;
end

% assessing if we want preconditioning
if strcmp(L,'off') || strcmp(L,'identity')
    precond = false;
    Lproject = false;
elseif ismatrix(L) && ~ischar(L)
    if size(L,2) ~= n
        error('The number of columns of the regularization matrix should be the same as the length of x')
    else
        if size(L,1) == n
            precond = true;
            Lproject = false;
        else
            precond = false;
            Lproject = true;
            p = size(L, 1);
            Lk = zeros(p,max(K));
        end
    end
else
    precond = true;
    Lproject = false;
    if strcmpi(L, 'Laplacian1D')
        L = LaplacianMatrix1D(n);
    elseif strcmpi(L, 'Laplacian2D')
        L = LaplacianMatrix2D(n);
    else
        error('Invalid string for regularization matrix')
    end
end

% Declare matrices.
X                = zeros(n,length(K));
Xnrm             = zeros(max(K),1);
Rnrm             = zeros(max(K),1);
RegParamVect     = zeros(max(K),1);
B                = zeros(max(K)+1,max(K)); 
V                = zeros(n, max(K));
rhs              = zeros(max(K)+1,1); % projected right-hand side
if restart
    saved_iterations = zeros(1, length(Ktot));
else
    saved_iterations = zeros(1, length(K));
end
GCV   = zeros(max(K),1);
Omega = zeros(max(K), 1);
warningGCV          = 0;
if notrue
    errornorms = false;
else
    errornorms = true;
    Enrm       = zeros(max(K),1);
    nrmtrue = norm(x_true(:));
    BestReg.RegP = [];
    BestReg.It = [];
    BestReg.X =[];
    BestReg.Enrm = [];
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestEnrm = 1e10;
end
% Main Code Begins Here
u = r;
U(:,1) = u/beta; % no matter the preconditioning (because it's on the right)
rhs(1) = beta;
% Iterate
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
    w = Atransp_times_vec(A, U(:,k));
    if precond
        % w = L' \ w;
        w = Ptransp_solve(L, w);
    end
    if k>1
        w = w - beta*V(:,k-1);
    end
    if reorth
        for jj = 1:k-1
            w = w - (V(:,jj)'*w)*V(:,jj);
        end
    end
    alpha = norm(w);
    V(:,k) = w/alpha;
    u = V(:,k);
    if precond
        % u = L \ u;
        u = P_solve(L, u);
    end
    u = A_times_vec(A, u);
    u = u - alpha*U(:,k);
    if reorth
        for jj = 1:k-1
            u = u - (U(:,jj)'*u)*U(:,jj);
        end
    end
    beta = norm(u);
    U(:,k+1) = u/beta;
    B(k,k) = alpha;
    B(k+1,k) = beta;
    rhsk = rhs(1:k+1); % current projected rhs
    
    
    if abs(alpha) <= eps || abs(beta) <= eps
        if verbose
            disp('Golub-Kahan bidiagonalization breaks down')
        end
        B = B(1:k+1,1:k);
        V = V(:,1:k);
        U = U(:,1:k+1);
        X(:,j+1) = x;
        X = X(:,1:j+1);
        if restart
            saved_iterations(j+1) = ktotcount-1;
        else
            saved_iterations(j+1) = k-1;
        end
        saved_iterations = saved_iterations(1:j+1);
        if k>1
            Xnrm    = Xnrm(1:k-1);
            Rnrm    = Rnrm(1:k-1);
            RegParamVect    = RegParamVect(1:k-1);
            if errornorms, Enrm = Enrm(1:k-1); end
        end
        % stop because the bidiagonal matrix is (numerically) singular
        % No chioce: even if NoStop is 'on'...we simpy cannot compute the solution, anymore
        if StopIt == MaxIter
            StopFlag = 'Breakdown of the Golub-Kahan bidiagonalization algorithm';
            StopReg.X = x; 
            StopReg.It = k-1;
            StopReg.RegP = RegParamk;
            StopReg.Xnrm = Xnrm(k-1);
            StopReg.Rnrm = Rnrm(k-1);
            if errornorms, StopReg.Enrm = Enrm(k-1); end
        end
        break
    end
    Bk = B(1:k+1,1:k);
    [Uk, Sk, Vk] = svd(Bk);
    if k==1
        Sk = Sk(1,1);
    else
        Sk = diag(Sk);
    end
    rhskhat = Uk'*rhsk;
    lsqr_res = abs(rhskhat(k+1))/nrmb;
    if Lproject
        Lk(:,k) = L*V(:,k);
        % update the Householder-QR factorization of Lk
        if k == 1
            [LUk, LRk] = householderQR(Lk(:,1:k));
        else
            [LUk, LRk] = upd_householderQR(Lk(:,1:k-1),...
            Lk(:,k), LUk, LRk);
        end
        LRksq = LRk(1:k,1:k);
        [Uk, Vk, ~, Ck, Sk] = gsvd(Bk, LRksq);
        rhskhat = Uk'*rhsk;
        if k==1
            gammak = Ck(1)/Sk(1);
        else
            gammak = sqrt(diag(Ck'*Ck)./diag(Sk'*Sk));
        end
    end
    
    if tik
        if isscalar(RegParam)
            RegParamk = RegParam;
            RegParamVect(k) = RegParamk;
        elseif strcmp(RegParam,'discrep')
            if k==1 
                RegParamVect(k) = RegParamk;
            end
        elseif strcmp(RegParam,'wgcv')
            if k>1
            if ~Lproject
                if adaptWGCV 
                    %Use the adaptive, weighted GCV method
                    Omega(k) = min(1, findomega(rhskhat, Sk));
                    omega = mean(Omega(1:k));
                end
                RegParamk = fminbnd('TikGCV', 0, Sk(1), [], rhskhat, Sk, omega);
                GCValk = GCVstopfun(RegParamk, Uk(1,:)', Sk, nrmb, m, n);
                
            else
                if adaptWGCV 
                    %Use the adaptive, weighted GCV method
                    Omega(k) = min(1, findomega(rhskhat, gammak));
                    omega = mean(Omega(1:k));
                end
                RegParamk = fminbnd('TikGCV', 0, gammak(k), [], rhskhat, gammak, omega);
                GCValk = GCVstopfun(RegParamk, Uk(1,:)', gammak, nrmb, m, n);
            end
            RegParamVect(k) = RegParamk;
            GCV(k) = GCValk;
            else
            RegParamk = 0; GCValk = 0;
            RegParamVect(k) = RegParamk; GCV(k) = GCValk;
            end
        elseif strcmp(RegParam,'gcv')
            if ~Lproject
                RegParamk = fminbnd('TikGCV', 0, Sk(1), [], rhskhat, Sk);
                GCValk = GCVstopfun(RegParamk, Uk(1,:)', Sk, nrmb, m, n);
            else
                RegParamk = fminbnd('TikGCV', 0, gammak(k), [], rhskhat, gammak);
                GCValk = GCVstopfun(RegParamk, Uk(1,:)', gammak, nrmb, m, n);
            end
            RegParamVect(k) = RegParamk;
            GCV(k) = GCValk;
        elseif strcmp(RegParam,'modgcv')
            if ~Lproject
                RegParamk = fminbnd('TikGCV', 0, Sk(1), [], rhskhat, Sk, m);
                GCValk = GCVstopfun(RegParamk, Uk(1,:)', Sk, nrmb, m, n);
            else
                RegParamk = fminbnd('TikGCV', 0, gammak(k), [], rhskhat, gammak, m);
                GCValk = GCVstopfun(RegParamk, Uk(1,:)', gammak, nrmb, m, n);
            end
            RegParamVect(k) = RegParamk;
            GCV(k) = GCValk;
        else
            error('Invalid parameter choice method')
        end
        if ~Lproject
            Dk = Sk.^2 + RegParamk^2;
            rhskhat = Sk .* rhskhat(1:k);
            yhat = rhskhat(1:k)./Dk;
            y = Vk * yhat;
        else
            BLk = [Bk; RegParamk*LRksq(1:k,:)];
            rhsLk = [rhsk; zeros(k,1)];
            y = BLk\rhsLk;
        end
        Rnrm(k) = norm(rhsk - Bk*y)/nrmb;
        d = V(:,1:k)*y;
        if precond, d = P_solve(L, d); end
        x = x0 + d;
        % Compute norms
        Xnrm(k) = norm(x(:));
        if errornorms
            Enrm(k) = norm(x_true(:) - x(:))/nrmtrue;
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
        AlreadySaved = 0;
        if any(k==K)
            j = j+1;
            X(:,j) = x;  
            saved_iterations(j) = k;
            % this is used to save the last iteration, in the case 
            % K = MaxIterIn
            % (when performing restarts, and the inner stopping criterion is not satisfied)
            if restart, saved_iterations(j) = ktotcount; end
            AlreadySaved = 1;              
        end
        if restart % else %
            if any(ktotcount == Ktot) && ~ AlreadySaved
                j = j+1;
                X(:,j) = x;
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
                StopReg.X = x;
                StopReg.It = k;
                StopReg.RegP = RegParamk;  
                StopReg.Xnrm = Xnrm(k);
                StopReg.Rnrm = Rnrm(k);
                if errornorms
                    Enrm = Enrm(1:k);
                    StopReg.Enrm = Enrm(k);
                end
                Xnrm    = Xnrm(1:k);
                Rnrm    = Rnrm(1:k);
                RegParamVect    = RegParamVect(1:k);
                V = V(:,1:k);
                B = B(1:k+1,1:k);
                U = U(:,1:k+1);
                X = X(:,1:j);
                saved_iterations = saved_iterations(1:j);
                if verbose
                    disp('reached maximum number of iterations')
                end
                StopFlag = 'reached maximum number of iterations';
                break
            end
        end       
        % update parameters, check stopping criteria
        if strcmp(RegParam,'discrep')
            if Rnrm(k) < eta*NoiseLevel
                % stopping criterion
                if StopIt == MaxIter % the method has not stopped, yet
                    if verbose
                        disp('The discrepancy principle is satisfied')
                    end
                    StopFlag = 'discrepancy principle (secant update method) satisfied';
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
                    StopReg.X = x;
                    StopReg.It = k;
                    StopReg.RegP = RegParamk;
                    StopReg.Xnrm = Xnrm(k);
                    StopReg.Rnrm = Rnrm(k);
                    if errornorms, StopReg.Enrm = Enrm(k); end
                    if ~ NoStop
                        Xnrm    = Xnrm(1:k);
                        Rnrm    = Rnrm(1:k);
                        RegParamVect    = RegParamVect(1:k);
                        V = V(:,1:k);
                        B = B(1:k+1,1:k);
                        U = U(:,1:k+1);
                        if errornorms, Enrm = Enrm(1:k); end
                        X = X(:,1:j);
                        saved_iterations = saved_iterations(1:j);
                        % stop because the discrepancy principle is satisfied
                        break
                    else
                        RegParamk = abs((eta*NoiseLevel - lsqr_res)/(Rnrm(k) - lsqr_res))*(RegParamk^2);
                        RegParamk = sqrt(RegParamk);
                        if k~=MaxIter, RegParamVect(k+1) = RegParamk; end
                    end
                else
                    RegParamk = abs((eta*NoiseLevel - lsqr_res)/(Rnrm(k) - lsqr_res))*(RegParamk^2);
                    RegParamk = sqrt(RegParamk);
                    if k~=MaxIter, RegParamVect(k+1) = RegParamk; end
                end
            else
                RegParamk = abs((eta*NoiseLevel - lsqr_res)/(Rnrm(k) - lsqr_res))*(RegParamk^2);
                RegParamk = sqrt(RegParamk);
                if k~=MaxIter, RegParamVect(k+1) = RegParamk; end
            end
        elseif strcmp(RegParam,'wgcv') || strcmp(RegParam,'gcv') || strcmp(RegParam,'modgcv')
            % check the stopping criterion (all the possibilities)
            if k > 1
            if strcmpi(stopGCV, 'GCVvalues')
                if abs((GCV(k)-GCV(k-1)))/GCV(2) < degflat && StopIt == MaxIter
                % the method has not stopped, yet
                        if verbose
                            disp('The stopping criterion for GCV principle is satisfied')
                        end
                        % stop because the GCV curve is too flat
                        StopFlag = 'GCV curve too flat';
                        if ~AlreadySaved && ~ NoStop
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
                        StopReg.X = x;
                        StopReg.It = k;
                        StopReg.RegP = RegParamk;
                        StopReg.Xnrm = Xnrm(k);
                        StopReg.Rnrm = Rnrm(k);
                        if errornorms, StopReg.Enrm = Enrm(k); end
                        if ~ NoStop
                            Xnrm    = Xnrm(1:k);
                            Rnrm    = Rnrm(1:k);
                            RegParamVect    = RegParamVect(1:k);
                            V = V(:,1:k);
                            B = B(1:k+1,1:k);
                            U = U(:,1:k+1);
                            if errornorms, Enrm = Enrm(1:k); end
                            X = X(:,1:j);
                            saved_iterations = saved_iterations(1:j);
                            % stop because the GCV stopping criterion is satisfied
                            break
                        end
                elseif GCV(k-1) < GCV(k) && ~ warningGCV && StopIt == MaxIter % Potential minimum reached. 
                    warningGCV = 1;
                    % Save data just in case.
                    x_save = x;
                    k_save = k; % for computing the GCV stopping criterion
                    j_save = j;
                    AlreadySaved_save = AlreadySaved;
                    RegParamk_save = RegParamk;
                    if restart, ktotcount_save = ktotcount; end
                elseif warningGCV && k > min(k_save + mintol, MaxIter) && StopIt == MaxIter % Passed window
                    if GCV(k_save) < GCV(k_save+1:min(k_save + mintol, MaxIter))
                        if verbose
                            disp('The stopping criterion for GCV principle is satisfied')
                        end
                        StopFlag = 'increasing GCV minima';
                        StopIt = k_save;
                        StopReg.It = k_save;
                        StopReg.X = x_save;
                        StopReg.RegP = RegParamk_save;
                        StopReg.Xnrm = Xnrm(k_save);
                        StopReg.Rnrm = Rnrm(k_save);
                        if errornorms
                            StopReg.Enrm = Enrm(k_save);
                        end
                        if ~ NoStop
                            j = j_save;
                            saved_iterations = saved_iterations(1:j);
                            X = X(:,1:j);
                            if ~AlreadySaved_save
                                j = j+1;
                                X(:,j) = x_save;
                                if restart
                                    saved_iterations(j) = ktotcount_save;
                                else
                                    saved_iterations(j) = k_save;
                                end
                            end
                            Xnrm    = Xnrm(1:k_save);
                            Rnrm    = Rnrm(1:k_save);
                            RegParamVect    = RegParamVect(1:k_save);
                            V = V(:,1:k_save);
                            B = B(1:k_save+1,1:k_save);
                            U = U(:,1:k_save+1);
                            if errornorms
                                Enrm = Enrm(1:k_save);
                                if BestReg.It > k_save
                                    [BestReg.Enrm, BestReg.It] = min(Enrm);
                                    BestReg.RegP = RegParamVect(BestReg.It);
                                    BestReg.Xnrm = Xnrm(BestReg.It);
                                    BestReg.Rnrm = Rnrm(BestReg.It);
                                    % recompute the best solution again
                                    ktemp = BestReg.It;
                                    Btemp = B(1:ktemp+1,1:ktemp);
                                    Vtemp = V(:,1:ktemp);
                                    rhsk = rhs(1:ktemp+1);
                                    [Uk, Sk, Vk] = svd(Btemp);
                                    if ktemp==1
                                        Sk = Sk(1,1);
                                    else
                                        Sk = diag(Sk);
                                    end
                                    rhskhat = Uk'*rhsk;
                                    if Lproject
                                        Lk = L*Vtemp;
                                        [~, LRksq] = qr(Lk,0);
                                    end
                                    if ~Lproject
                                        Dk = Sk.^2 + RegParamk^2;
                                        rhskhat = Sk .* rhskhat(1:ktemp);
                                        yhat = rhskhat(1:ktemp)./Dk;
                                        y = Vk * yhat;
                                    else
                                        HLk = [Htemp; RegParamk*LRksq];
                                        rhsLk = [rhsk; zeros(ktemp,1)];
                                        y = HLk\rhsLk;
                                    end
                                    dtemp = Vtemp*y;
                                    if precond, dtemp = P_solve(L, dtemp); end
                                    xtemp = x0 + dtemp;
                                    BestReg.X = xtemp;
                                end
                            end
                            X = X(:,1:j);
                            saved_iterations = saved_iterations(1:j);
                            if restart, ktotcount = ktotcount_save; end
                            % stop because the GCV stopping criterion is satisfied
                            break
                        end
                    else
                        warningGCV = 0;
                    end
                end
            elseif strcmpi(stopGCV, 'resflat')
                if abs((Rnrm(k)-Rnrm(k-1)))/Rnrm(k-1) < resdegflat && ...
                    Rnrm(k) == min(Rnrm(1:k)) && StopIt == MaxIter
                    if verbose
                        disp('The stopping criterion for GCV principle is satisfied')
                    end
                    % stop because discrepancy (i.e., residual for the
                    % regularized problem) stabilizes
                    StopFlag = 'the residual norm stabilizes';
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
                    StopReg.X = x;
                    StopReg.It = k;
                    StopReg.RegP = RegParamk;
                    StopReg.Xnrm = Xnrm(k);
                    StopReg.Rnrm = Rnrm(k);
                    if errornorms, StopReg.Enrm = Enrm(k); end
                    if ~ NoStop
                        Xnrm    = Xnrm(1:k);
                        Rnrm    = Rnrm(1:k);
                        RegParamVect    = RegParamVect(1:k);
                        V = V(:,1:k);
                        B = B(1:k+1,1:k);
                        U = U(:,1:k+1);
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
end
if k == MaxIter 
    if StopIt == MaxIter
        % Stop because max number of iterations reached
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
        StopReg.X = x;
        StopReg.It = k;
        StopReg.RegP = RegParamk;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        if errornorms, StopReg.Enrm = Enrm(k); end
        Xnrm    = Xnrm(1:k);
        Rnrm    = Rnrm(1:k);
        RegParamVect    = RegParamVect(1:k);
        V = V(:,1:k);
        B = B(1:k+1,1:k);
        U = U(:,1:k+1);
        if errornorms, Enrm = Enrm(1:k); end
        X = X(:,1:j);
        saved_iterations = saved_iterations(1:j);
    end 
end
if ~noIterBar, close(h_wait), end
if nargout==2
  if NoStop
      info.its = k;
  else
      info.its = StopIt;
  end
  info.saved_iterations = saved_iterations(1:j);
  info.StopReg = StopReg;
  info.StopFlag = StopFlag;
  if errornorms
    info.Enrm = Enrm;
    info.BestReg = BestReg;
  end
  info.Xnrm = Xnrm;
  info.Rnrm = Rnrm;
  info.RegP = RegParamVect;
  if strcmp(RegParam,'wgcv') || strcmp(RegParam,'gcv') || strcmp(RegParam,'modgcv')
    info.GCValues = GCV(1:k);
  end
  if strcmp(DecompOut,'on')
      info.V = V;
      info.U = U;
      info.B = B;
  end
  if restart
      info.ktotcount = ktotcount;
  end
end

%% ---------------SUBFUNCTIONS ---------------------------------------
function [U,R] = householderQR(L)
%   
%  [U,R] = householderQR(L)
%  This function computes the Householder-QR factorization of L
%  (a "projected" regularization matrix), that will be used to define a
%  "projected" regularization matrix R to employ within the LSQR iterates.
%

[m, n] = size(L);
R = L;
U = zeros(m, n);
for k = 1:n
    x = L(k:m,k);
    e = zeros(length(x),1); e(1) = 1;
    u = sign(x(1))*norm(x(:))*e + x;
    u = u./norm(u(:));
    R(k:m, k:n) = R(k:m, k:n) -2*u*(u'*R(k:m, k:n));
    U(k:m,k) = u;
end

function [U,R] = upd_householderQR(L,ll,U,R)
%   
% [U,R] = upd_householderQR(L, ll, U, R)
% This function updates the Householder-QR factorization of [L, ll].
%
% Input:
%   L  - matrix whose QR factorization is defined by U and R
%   ll - column appended to L, i.e., [L, ll]
%    U - matrix defining the orthogonal matrix Q, such that L = QR
%    R - upper triangular factor of L = QR
%

[m,n] = size(L);
Unew = zeros(m, n+1);
Unew(:,1:n) = U;
w = ll;
for i = 1:n
    u = U(i:m,i);
    w(i:m) = w(i:m) - 2*u*(u'*w(i:m));
end
v = w(1:n); x = w(n+1:m);
e = zeros(length(x),1); e(1) = 1;
u = sign(x(1))*norm(x(:))*e + x;
u = u./norm(u(:));
x = x -2*u*(u'*x);
Unew(n+1:m,n+1)=u;
U = Unew;
rr = [v; x];
R = [R, rr];

function omega = findomega(bhat, s)
%
%   omega = findomega(bhat, s)
%
%  This function computes a value for the omega parameter used in wGCV.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the 
%  derivative equal to zero and then solve for omega.
%  
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = length(s);

t0 = sum(abs(bhat(n+1:m)).^2);
alpha = s(end);
s2 = abs(s) .^ 2;
alpha2 = alpha^2;

tt = 1 ./ (s2 + alpha2);

t1 = sum(s2 .* tt);
t2 = abs(bhat(1:n).*alpha.*s) .^2;
t3 = sum(t2 .* abs((tt.^3)));

t4 = sum((s.*tt) .^2);
t5 = sum((abs(alpha2*bhat(1:n).*tt)).^2);

v1 = abs(bhat(1:n).*s).^2;
v2 = sum(v1.* abs((tt.^3)));

%
% Now compute omega.
%
omega = (m*alpha2*v2)/(t1*t3 + t4*(t5 + t0));
