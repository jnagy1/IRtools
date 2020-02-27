function [X,info] = IRcgls(A,b,varargin)
% IRcgls Conjugate Gradient algorithm for Least Squares problems
%
% options  = IRcgls('defaults')
% [X,info] = IRcgls(A,b)
% [X,info] = IRcgls(A,b,K)
% [X,info] = IRcgls(A,b,options)
% [X,info] = IRcgls(A,b,K,options)
%
% This function applies the CG algorithm implicitly to the normal equations
% for the least squares problem.  We obtain a regularized solution by
% terminating the iterations.  Alternatively this function can be used to
% computed a Tikhonov solution.
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
%      MaxIter    - maximum allowed number of iterations
%                   [ positive integer | {100} ]
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
%      RegParam   - regularization parameter lambda, to be employed if
%                   CGLS is used to solve the regularized problem
%                   (A'*A + lambda^2*L'*L)*x = A'*b;
%                   [ {0} | nonnegative scalar ]
%      RegMatrix  - regularization matrix L, used either as a
%                   regularization matrix to solve the regularized
%                   problem (A'*A + lambda^2*L'*L)*x = A'*b, or as a
%                   priorconditioner for the problem (in this case
%                   one should set RegParam = 0)
%                   [ {'Identity'} | 'Laplacian1D' | 'Laplacian2D' |
%                     matrix | function handle ]
%      Reorth     - indicates if reorthogonalization should be applied
%                   [ {'on'} | 'off' ]
%      IterBar    - shows the progress of the iterations
%                   [ {'on'} | 'off' ]
%      NoStop     - specifies whether the iterations should proceed after
%                   a stopping criterion has been satisfied
%                   [ 'on' | {'off'} ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its      - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - a string that describes the stopping condition:
%                   * Reached maximum number of iterations
%                   * Residual tolerance satisfied (discrepancy principle) 
%                   * Normal equation residual tolerance satisfied
%      StopReg  - struct containing information about the solution that
%                 satisfies the stopping criterion, with the fields:
%                   It   : iteration where the stopping criterion is satisfied
%                   X    : the solution satisfying the stopping criterion
%                   Enrm : the corresponding relative error (requires x_true)
%      Rnrm     - relative residual norms at each iteration
%      NE_Rnrm  - normal eqs relative residual norms
%      Xnrm     - solution norms at each iteration
%      Enrm     - relative error norms (requires x_true) at each iteration
%      BestReg  - struct containing information about the solution that
%                 minimizes Enrm (requires x_true), with the fields:
%                   It   : iteration where the minimum is attained
%                   X    : best solution
%                   Enrm : best relative error
%      StdCGLS  - struct containing information about the standard CGLS,
%                 in case enriched CGLS is used
%
% See also: IRenrich, IRhybrid_lsqr, IRnnfcgls, IRrrgmres, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0','none', 'MaxIter',100 , 'x_true','none',...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, 'IterBar','on', ...
    'stdCGLS_out','off', 'NoStop','off', ...
    'RegParam',0, 'RegMatrix','Identity', 'Reorth','on');
  
% If input is 'defaults,' return the default options in X.
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    X = defaultopt;
    return;
end

defaultopt.restart    = 'off';
defaultopt.verbosity  = 'on';
defaultopt.enrichment = 'none';
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
        if isfield(options, 'MaxIter') && ~isempty(options.MaxIter) && ... 
                (~isempty(K) && options.MaxIter ~= max(K))
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
Reorth     = IRget(options, 'Reorth',     [], 'fast');
L          = IRget(options, 'RegMatrix',  [], 'fast');
TikParam   = IRget(options, 'RegParam',   [], 'fast');
NoStop     = IRget(options, 'NoStop',     [], 'fast');
Q_tilde    = IRget(options, 'enrichment', [], 'fast');
restart    = IRget(options, 'restart',    [], 'fast');
verbose    = IRget(options, 'verbosity',  [], 'fast');

restart = strcmp(restart, 'on');
verbose = strcmp(verbose, 'on');

if ischar(TikParam), TikParam = 0; end

if isempty(K)
    K = MaxIter;
end
% Sorting the iteration numbers (in case they are shuffled in input).
K = K(:); K = sort(K,'ascend'); K = unique(K);
if ~((isreal(K) && (all(K > 0)) && all(K == floor(K))))
    error('K must be a vector of positive real integers')
end
if K(end) ~= MaxIter
    MaxIter = K(end); 
end

if strcmp(Reorth,'on')
    reorth = true; 
else
    reorth = false;
end

StopIt = MaxIter;

%  We need to find the number of columns in matrix A, but if A is not given 
%  as a matrix, and no initial guess is given, then we can find it by 
%  computing A'*b.  We need this anyway, so doesn't cost additional work.

d = Atransp_times_vec(A, b);
n = length(d);
m = length(b);

nrmb = norm(b(:));
nrmAtb = norm(d(:));

if restart
    ktotcount  = IRget(options, 'ktotcount', [],'fast');
    TotIterMax = IRget(options, 'TotIterMax',[],'fast');
    if strcmp(TotIterMax, 'none') || TotIterMax < MaxIter
        TotIterMax = MaxIter;
    end
    if strcmp(ktotcount, 'none')
        error('The total iteration counter must be assigned')
    end
    Ktot = IRget(options, 'Ktot', [], 'fast');
    % No checks on Ktot, it should be given from IRrestart.
end

if isempty(NoiseLevel) || strcmp(NoiseLevel,'none')
    Rtol = 0;
else
    Rtol = eta*NoiseLevel;
end

% See if an initial guess is given. If not, then use 0 as the initial guess.  
x = IRget(options, 'x0', [], 'fast');

if strcmp(x,'none')
    r = b;
    x = zeros(n,1);
else
    try
        Ax = A_times_vec(A, x);
        r = b(:) - Ax;
        d = d - Atransp_times_vec(A, Ax);
    catch
        error('Check the length of x')
    end
end

% Enrichment?
if ischar(Q_tilde)
    if strcmpi(Q_tilde,'none')
        enrich = false;
    else
        enrich = true;
        Q_tilde = ones(n,1)/sqrt(n);
        warning(['The default enrichment subspace is spanned by a constant vector. ',...
                 'Note that this may not be a good choice for your specific problem.'])
    end
else
    enrich = true;
end

% Tikhonov regularization?
% If we perform enrichment it is possible to use standard-form regulariza-
% tion; this is done in a different manner than for CGLS.  Otherwise:
%
% If TikParam is set to 'off', then any regularization operators are 
% implemented through preconditioning.  In this case, a regularization
% parameter is not specified, and regularization is enforced through 
% termination of the iteration (i.e., we use iterative regularization).
%
% If TikParam is 'on', then the least squares system is augmented
% with the regularization operator, weighted by the given regularization
% parameter, and the resulting over-determined least squares problem is
% solved.  This is the usual implementation of Tikhonov regularization.
enrichTik = false;
if enrich && TikParam ~= 0
    enrichTik = true;
    tik = false;
    precond = false;
elseif TikParam == 0
    tik = false;
    if strcmp(L,'off') || strcmp(L,'identity')
        precond = false;
    else
        precond = true;
        if strcmpi(L, 'Laplacian1D')
            L = LaplacianMatrix1D(n);
        elseif strcmpi(L, 'Laplacian2D')
            L = LaplacianMatrix2D(n);
        %else
        % Assume the user has supplied a matrix of function handle for L.
        end
    end
else
    precond = false;
    tik = true;
    if strcmpi(L,'identity')
        L = speye(n);
    elseif strcmpi(L, 'Laplacian1D')
        L = LaplacianMatrix1D(n);
    elseif strcmpi(L, 'Laplacian2D')
        L = LaplacianMatrix2D(n);
    %else
    %   Assume a user has supplied a matrix of function handle for L.
    end
    Lx = A_times_vec(L, TikParam*x);
    r = [ r ; -Lx ];
    d = d - Atransp_times_vec(L, TikParam*Lx);
end

if (enrichTik || tik) && Rtol ~= 0
    warning('With these input options IRcgls solves the normal equations associated to the Tikhonov regularized problem. The solution to this problem should be computed with high accuracy, so ''NoiseLevel'' should ideally be ''none'' or 0.')
end

% Declare matrices.
X = zeros(n,length(K));
if reorth, Q = zeros(n,MaxIter+1); end
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
    BestReg.Xnrm = [];
    BestReg.Rnrm = [];
    BestReg.NE_Rnrm = [];
    BestReg.Enrm = [];
    BestEnrm = 1e10;
end

if enrich
    [q1, q2] = size(Q_tilde);
    if q1 ~= n
        error('The columns of the enrichment matrix must have the same length as the solution vector')
    end
    V = zeros(m, q2);
    S = zeros(n, q2);
    for i = 1:q2
        V(:,i) = A_times_vec(A, Q_tilde(:,i));
        S(:,i) = Atransp_times_vec(A, V(:,i));
    end
    stdCGLS_out = IRget(options, 'stdCGLS_out', [], 'fast');
    stdCGLS_out = strcmp(stdCGLS_out, 'on');
    X_tilde       = zeros(n,length(K));
    Xnrm_tilde    = zeros(MaxIter,1);
    Rnrm_tilde    = zeros(MaxIter,1);
    NE_Rnrm_tilde = zeros(MaxIter,1);
    if errornorms
        Enrm_tilde = zeros(MaxIter,1);
        BestReg_tilde.It = [];
        BestReg_tilde.X = [];
        BestReg_tilde.Xnrm = [];
        BestReg_tilde.Rnrm = [];
        BestReg_tilde.NE_Rnrm = [];
        BestReg_tilde.Enrm = [];
        BestEnrm_tilde = 1e10;
    end
end

NoStop = strcmp(NoStop,'on');

if restart
    saved_iterations = zeros(1, length(Ktot));
else
    saved_iterations = zeros(1, length(K));
end

% Prepare for iterations.
if precond
    d = Ptransp_solve(L, d);
end
if reorth, Q(:,1) = d/norm(d(:)); end
normr2 = norm(d(:))^2;
if precond
    d = P_solve(L, d);
end
if enrichTik
    s = Atransp_times_vec(A, r) - TikParam^2*x;
end

% Iterate.
noIterBar = strcmp(IterBar,{'off'});
if ~noIterBar
    h_wait = waitbar(0, 'Running iterations, please wait ...');
end
j = 0;
for k=1:MaxIter
    if restart, ktotcount = ktotcount + 1; end
    if ~noIterBar
        waitbar(k/MaxIter, h_wait)
    end
    
    % Update x and r vectors.
    if tik
        Ad = [A_times_vec(A, d); A_times_vec(L, TikParam*d)];
    else
        Ad = A_times_vec(A, d);
    end
    if enrich
        f = Atransp_times_vec(A, Ad);
    end
    normAd2 = Ad'*Ad;
    if enrichTik
        alpha = normr2/(normAd2 + TikParam^2*(d'*d));
    else
        alpha = normr2/normAd2;
    end
    x  = x + alpha*d;
    r  = r - alpha*Ad;
    if tik
        s = Atransp_times_vec(A, r(1:m)) + Atransp_times_vec(L, TikParam*r(m+1:end));
    elseif enrichTik
        s = s - alpha*(f + TikParam^2*d);
    else
        s  = Atransp_times_vec(A, r);
    end
    if precond
        q = Ptransp_solve(L, s);
    else
        q = s;  
    end
    if reorth
        for i=1:k, q = q - (Q(:,i)'*q)*Q(:,i); end
        Q(:,k+1) = q/norm(q);
    end
    % normr2_new = norm(s)^2; %%% there was a bug!!!
    normr2_new = norm(q)^2;
    if precond
        s = P_solve(L, q);
    else
        s = q;    
    end
        
    % Update d vector.
    beta = normr2_new/normr2;
    normr2 = normr2_new;
    if enrich
        d_old = d;
    end
    d = s + beta*d;
    
    if enrich
        if enrichTik
            g = (V'*Ad + TikParam^2*(Q_tilde'*d_old))/...
                (normAd2 + TikParam^2*(d_old'*d_old));
        else
           g = (V'*Ad)/normAd2;
        end
        Q_tilde  = Q_tilde - d_old*g';
        V = V - Ad*g';
        % Ought to update QR factorization, but \ is faster in Matlab.
        if enrichTik
            y_tilde = [V;TikParam*Q_tilde]\[zeros(m,1);s/TikParam];
        else
            y_tilde = V\r;
        end
        x_tilde = x + Q_tilde*y_tilde;
        r_tilde = r - V*y_tilde;
        S = S - f*g';
        if enrichTik
            d_tilde = d - (S + TikParam^2*Q_tilde)*y_tilde;
        else
            d_tilde = d - S*y_tilde;
        end
    end
    
    % Compute norms.
    Xnrm(k)    = norm(x(:));
    Rnrm(k)    = norm(r(:))/nrmb;
    NE_Rnrm(k) = norm(s(:))/nrmAtb;
    if enrich
        Xnrm_tilde(k)    = norm(x_tilde(:));
        Rnrm_tilde(k)    = norm(r_tilde(:))/nrmb;
        NE_Rnrm_tilde(k) = norm(d_tilde(:))/nrmAtb;
    end
    if errornorms
        Enrm(k) = norm(x_true(:)-x(:))/nrmtrue;
        if Enrm(k)<BestEnrm
            BestReg.It = k;
            BestReg.X = x;
            BestReg.Xnrm = Xnrm(k);
            BestReg.Rnrm = Rnrm(k);
            BestReg.NE_Rnrm = NE_Rnrm(k);
            BestEnrm = Enrm(k);
            BestReg.Enrm = BestEnrm;
        end
        if enrich
            Enrm_tilde(k) = norm(x_true(:)-x_tilde(:))/nrmtrue;
            if Enrm_tilde(k)<BestEnrm_tilde
                BestReg_tilde.It = k;
                BestReg_tilde.X = x_tilde;
                BestReg_tilde.Xnrm = Xnrm_tilde(k);
                BestReg_tilde.Rnrm = Rnrm_tilde(k);
                BestReg_tilde.NE_Rnrm = NE_Rnrm_tilde(k);
                BestEnrm_tilde = Enrm_tilde(k);
                BestReg_tilde.Enrm = BestEnrm_tilde;
            end
        end
    end
    AlreadySaved = 0;
    if any(k==K)
        j = j+1;
        X(:,j) = x;
        if enrich
            X_tilde(:,j) = x_tilde;
        end
        saved_iterations(j) = k;
        % This is used to save the last iteration in the case K = MaxIterIn
        % (when performing restarts, and the inner stopping criterion is
        % not satisfied)
        if restart, saved_iterations(j) = ktotcount; end
        AlreadySaved = 1; 
    end
    if restart
        if any(ktotcount == Ktot) && ~ AlreadySaved
            j = j+1;
            X(:,j) = x;
            if enrich, X_tilde(:,j) = x_tilde; end
            saved_iterations(j) = ktotcount;
            AlreadySaved = 1;                
        end
        if ktotcount == TotIterMax
            if ~ AlreadySaved
                j = j+1;
                X(:,j) = x;
                if enrich, X_tilde(:,j) = x_tilde; end
                saved_iterations(j) = ktotcount;
                AlreadySaved = 1;
            end
            StopIt = k;
            StopReg.It = k;
            StopReg.X = x;
            StopReg.Xnrm = Xnrm(k);
            StopReg.Rnrm = Rnrm(k);
            StopReg.NE_Rnrm = NE_Rnrm(k);
            if enrich
                StopReg_tilde.It = k;
                StopReg_tilde.X = x_tilde;
                StopReg_tilde.Xnrm = Xnrm_tilde(k);
                StopReg_tilde.Rnrm = Rnrm_tilde(k);
                StopReg_tilde.NE_Rnrm = NE_Rnrm_tilde(k);
            end
            if errornorms
                Enrm = Enrm(1:k);
                if enrich, Enrm_tilde = Enrm_tilde(1:k); end
                StopReg.Enrm = Enrm(k);
                if enrich, StopReg_tilde.Enrm = Enrm_tilde(k); end
            end
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if enrich
                Xnrm_tilde    = Xnrm_tilde(1:k);
                Rnrm_tilde    = Rnrm_tilde(1:k);
                NE_Rnrm_tilde = NE_Rnrm_tilde(1:k);
            end
            if errornorms
                Enrm = Enrm(1:k);
                if enrich, Enrm_tilde = Enrm_tilde(1:k); end
            end
            X = X(:,1:j);
            if enrich, X_tilde = X_tilde(:,1:j); end
            saved_iterations = saved_iterations(1:j);
            if verbose
                disp('Reached maximum number of iterations')
            end
            StopFlag = 'Reached maximum number of iterations';
            break
        end
    end
    
    if enrich
    if (Rnrm_tilde(k) <= Rtol) && (StopIt == MaxIter)
        % Stop because residual satisfies ||b-A*x||/||b|| <= Rtol.
        if verbose
            disp('Enriched residual tolerance satisfied')
        end
        StopFlag = 'Enriched residual tolerance satisfied';
        if ~AlreadySaved && ~NoStop
            j = j+1;
            X(:,j) = x;
            X_tilde(:,j) = x_tilde;
            if restart
                saved_iterations(j) = ktotcount;
            else
                saved_iterations(j) = k;
            end
            AlreadySaved = 1;
        end
        StopIt = k;
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        StopReg_tilde.It = k;
        StopReg_tilde.X = x_tilde;
        StopReg_tilde.Xnrm = Xnrm_tilde(k);
        StopReg_tilde.Rnrm = Rnrm_tilde(k);
        StopReg_tilde.NE_Rnrm = NE_Rnrm_tilde(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
            StopReg_tilde.Enrm = Enrm_tilde(k);
        end
        if ~ NoStop
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            Xnrm_tilde    = Xnrm_tilde(1:k);
            Rnrm_tilde    = Rnrm_tilde(1:k);
            NE_Rnrm_tilde = NE_Rnrm_tilde(1:k);
            if errornorms
                Enrm       = Enrm(1:k);
                Enrm_tilde = Enrm_tilde(1:k);
            end
            X = X(:,1:j);
            X_tilde = X_tilde(:,1:j);
            saved_iterations = saved_iterations(1:j);
            break
        end
    end
    if (NE_Rnrm_tilde(k) <= NE_Rtol) && (StopIt == MaxIter)
        if verbose
            disp('Enriched normal equations residual tolerance satisfied')
        end
        StopFlag = 'Enriched normal equations residual tolerance satisfied';
        if ~AlreadySaved && ~NoStop
            j = j+1;
            X(:,j) = x;
            X_tilde(:,j) = x_tilde;
            if restart    
                saved_iterations(j) = ktotcount;
            else
                saved_iterations(j) = k;
            end
            AlreadySaved = 1;
        end
        StopIt = k;
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        StopReg_tilde.It = k;
        StopReg_tilde.X = x_tilde;
        StopReg_tilde.Xnrm = Xnrm_tilde(k);
        StopReg_tilde.Rnrm = Rnrm_tilde(k);
        StopReg_tilde.NE_Rnrm = NE_Rnrm_tilde(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
            StopReg_tilde.Enrm = Enrm_tilde(k);
        end
        if ~ NoStop
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            Xnrm_tilde    = Xnrm_tilde(1:k);
            Rnrm_tilde    = Rnrm_tilde(1:k);
            NE_Rnrm_tilde = NE_Rnrm_tilde(1:k);
            if errornorms
                Enrm       = Enrm(1:k);
                Enrm_tilde = Enrm_tilde(1:k);
            end
            X = X(:,1:j);
            X_tilde = X_tilde(:,1:j);
            saved_iterations = saved_iterations(1:j);
            break
        end       
    end 
    else
    if (Rnrm(k) <= Rtol) && (StopIt == MaxIter)
        % Stop because residual satisfies ||b-A*x||/||b|| <= Rtol.
        if verbose
            disp('Residual tolerance satisfied')
        end
        StopFlag = 'Residual tolerance satisfied';
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
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
        end
        if ~ NoStop
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if errornorms
                Enrm = Enrm(1:k);
            end
            X = X(:,1:j);
            if enrich, X_tilde = X_tilde(:,1:j); end
            saved_iterations = saved_iterations(1:j);
            break
        end
    end
    if (NE_Rnrm(k) <= NE_Rtol) && (StopIt == MaxIter)
        if verbose
            disp('Normal equations residual tolerance satisfied')
        end
        StopFlag = 'Normal equations residual tolerance satisfied';
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
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if errornorms
            StopReg.Enrm = Enrm(k);
        end
        if ~ NoStop
            Xnrm    = Xnrm(1:k);
            Rnrm    = Rnrm(1:k);
            NE_Rnrm = NE_Rnrm(1:k);
            if errornorms
                Enrm = Enrm(1:k);
            end
            X = X(:,1:j);
            if enrich, X_tilde = X_tilde(:,1:j); end
            saved_iterations = saved_iterations(1:j);
            break
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
        StopFlag = 'Reached maximum number of iterations';
        if ~AlreadySaved
            j = j+1;
            X(:,j) = x;
            if enrich, X_tilde(:,j) = x_tilde; end
            if restart
                saved_iterations(j) = ktotcount;
            else
                saved_iterations(j) = k;
            end
        end
        StopReg.It = k;
        StopReg.X = x;
        StopReg.Xnrm = Xnrm(k);
        StopReg.Rnrm = Rnrm(k);
        StopReg.NE_Rnrm = NE_Rnrm(k);
        if enrich
            StopReg_tilde.It = k;
            StopReg_tilde.X = x_tilde;
            StopReg_tilde.Xnrm = Xnrm_tilde(k);
            StopReg_tilde.Rnrm = Rnrm_tilde(k);
            StopReg_tilde.NE_Rnrm = NE_Rnrm_tilde(k);
        end
        if errornorms
            StopReg.Enrm = Enrm(k);
            if enrich, StopReg_tilde.Enrm = Enrm_tilde(k); end
        end
        Xnrm    = Xnrm(1:k);
        Rnrm    = Rnrm(1:k);
        NE_Rnrm = NE_Rnrm(1:k);
        if enrich
            Xnrm_tilde    = Xnrm_tilde(1:k);
            Rnrm_tilde    = Rnrm_tilde(1:k);
            NE_Rnrm_tilde = NE_Rnrm_tilde(1:k);
        end
        if errornorms
            Enrm = Enrm(1:k);
            if enrich, Enrm_tilde = Enrm_tilde(1:k); end
        end
        X = X(:,1:j);
        if enrich, X_tilde = X_tilde(:,1:j); end
        saved_iterations = saved_iterations(1:j);
    end 
end

if ~noIterBar, close(h_wait), end
if nargout==2
  info.its = k;
  if restart
      info.ktotcount = ktotcount;
  end
  info.saved_iterations = saved_iterations(1:j);
  if enrich
    StdCGLS.X = X;
    StdCGLS.StopReg = StopReg;
    StdCGLS.Xnrm = Xnrm;
    StdCGLS.Rnrm = Rnrm;
    StdCGLS.NE_Rnrm = NE_Rnrm;
    if errornorms
        StdCGLS.Enrm = Enrm;
        StdCGLS.BestReg = BestReg;
    end
    X       = X_tilde;
    StopReg = StopReg_tilde;
    Xnrm    = Xnrm_tilde(1:k);
    Rnrm    = Rnrm_tilde(1:k);
    NE_Rnrm = NE_Rnrm_tilde(1:k);
    if errornorms
        Enrm    = Enrm_tilde;
        BestReg = BestReg_tilde;
    end
  end
  info.StopFlag = StopFlag;
  info.StopReg = StopReg;
  info.Rnrm = Rnrm(1:k);
  info.NE_Rnrm = NE_Rnrm(1:k);
  info.Xnrm = Xnrm(1:k);
  if errornorms
    info.Enrm = Enrm(1:k);
    info.BestReg = BestReg;
  end
  if enrich && stdCGLS_out
      info.StdCGLS = StdCGLS;
  end
end