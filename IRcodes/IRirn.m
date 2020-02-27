function [X, info] = IRirn(A, b, varargin)
% IRirn Least squares solver with 1-norm penalization term
%
% options  = IRirn('defaults')
% [X,info] = IRirn(A,b)
% [X,info] = IRirn(A,b,K)
% [X,info] = IRirn(A,b,options)
% [X,info] = IRirn(A,b,K,options)
%
% Iteratively reweighted norm algorithm for computing a 1-norm penalized 
% solution. 
%
% IRirn is a simplified driver for IRrestart, which uses an inner-outer 
% iteration scheme. Semi-convergent or hybrid iterative solvers are used 
% in the inner iterations, using one of the iterative methods in IRtools 
% (e.g., IRhybrid_gmres). In the case of IRirn, the 1-norm penalization 
% is updated at each outer iteration. 
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
%  b : right-hand side vector
%  K : (optional) integer vector that specifies which (total) iterates are 
%      returned in X; the maximum number of iterations is assumed to be max(K)
%      [ positive integer | vector of positive components ]
%  options : structure with the following fields (optional)
%      x0            - initial guess for the iterations; default = zero vector
%                      [ array | {'none'} ]
%      MaxIterIn     - maximum number of inner iterations
%      MaxIterOut    - maximum number of outer iterations
%      x_true        - true solution; allows us to returns error norms with
%                      respect to x_true at each iteration
%                      [ array | {'none'} ]
%      RegParam      - a value or a method to find the regularization used
%                      in the inner iterations: 
%                      [ non-negative scalar | {'gcv'} | 'discrep' ]
%                      This also determines which stopping rule is used for
%                      the inner iterations.
%                      If 'gcv' is chosen, the inner iteration is stopped
%                        when the GCV function minimum stabilizes or increases 
%                        within a certain window of iterations (see 'stopGCV',
%                        'FlatTol' and 'MinTol').
%                      If 'discrep' is chosen, and NoiseLevel is provided,
%                        then the discrepancy principle is used as stopping
%                        criterion (see 'NoiseLevel' and 'eta').
%      stopGCV       - stopping criterion for the inner iterations when
%                      GCV is used
%                      [ GCVvalues | {'resflat'} ]
%      FlatTol       - tolerance for detecting flatness (stabilization)
%                      in the GCV function as a stopping criterion for the
%                      inner iterations
%                      [ {10^-6} | non-negative scalar ]
%      MinTol        - window of iterations: if the GCV minimum continues
%                      to increase over this window, then the inner
%                      iterations are stopped
%                      [ {3} | positive integer ]
%      RegMatrix     - priorconditioner for the inner iterations
%                      [ {'identity'} | square nonsingular matrix | 
%                        function handle ]
%      NoiseLevel    - norm of noise in rhs divided by norm of rhs (must be
%                      assigned if RegParam is 'discrep')
%                      [ {'none'} | nonnegative scalar ]
%      eta           - safety factor for the discrepancy principle
%                      [ {1.01} | scalar greater than (and close to) 1 ]
%      RegParam0     - first regularization parameter, used only on the 
%                      very first iteration (needed if RegParam is 'discrep')
%                      [ {1} | positive scalar ]
%      stopOut       - stopping criterion for the outer iterations;
%                      [ {'xstab'} | 'Lxstab' | 'regPstab' ]
%      inSolver      - solver to be employed during the inner iterations
%                      [ {'gmres'} | 'lsqr' | 'fgmres' | 'cgls']
%      adaptConstr   - approximate constraint or regularization to be
%                      incorporated
%                      [ {'sp'} | 'spnn' | 'none' ]
%      nonnegativity - may be used to also impose nonnegativity
%                      (similarly to 'spnn')
%                      [ 'on' | {'off'} ]
%      IterBar       - shows the progress of the outer iterations
%                      [ {'on'} | 'off' ]
%      NoStopIn      - specifies whether the inner iterations should
%                      proceed after a stopping criterion has been satisfied
%                      [ 'on' | {'off'} ]
%      NoStopOut     - specifies whether the outer iterations should
%                      proceed after a stopping criterion has been satisfied
%                      [ 'on' | {'off'} ]
%      verbosity     - switch on or off the "verbosity" of the function
%                      [ {'on'} | 'off' ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its          - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag_in  - string that describes the inner stopping condition:
%                       * Stopping criterion of the inner iterations is
%                         never satisfied
%                       * Stopping criterion is satisfied at least once
%                         during the inner iterations
%      StopFlag_out - string that describes the outer stopping condition;
%                     depending on the inputs it can be one of the following:
%                       * Outer stopping criterion is never satisfied
%                       * Diagonal weighting matrix is numerically zero
%                       * Solution stabilizes
%                       * Transformed solution stabilizes
%                       * Tegularization parameter stabilizes
%      Rnrm     - relative residual norms at each iteration
%      Xnrm     - solution norms at each iteration
%      Enrm     - relative error norms (requires x_true) at each iteration
%      Xout     - approximate solutions at the end of each inner cycle,
%                 stored column-wise
%      itsInOut - 3-column matrix whose the columns store
%                   1. outer iteration count
%                   2. inner iteration count (i.e., for each cycle)
%                   3. total iteration count
%      StopReg  - struct containing information about the solution that
%                 satisfies the stopping criterion.  Fields:
%                   It   : iteration where the stopping criterion is satisfied
%                   X    : solution satisfying the stopping criterion
%                   Enrm : the corresponding relative error (requires x_true)
%      BestReg  - struct containing information about the solution that
%                 minimizes Enrm (requires x_true). Fields:
%                   It   : iteration where the minimum is attained
%                   X    : best solution
%                   Enrm : best relative error
%
% See also: IRell1, IRhtv, IRhybrid_fgmres, IRrestart, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0', 'none', 'MaxIterIn', 30 , 'MaxIterOut', 20 , ...
    'RegParam', 'gcv', 'stopGCV', 'resflat', ...
    'resflatTol', 0.05, 'GCVflatTol', 10^-6, 'GCVminTol', 3,...
    'x_true', 'none', 'IterBar', 'on', 'NoStop', 'off', 'NoStopIn', 'off',...
    'NoStopOut', 'off', 'stopOut', 'xstab', 'stabOut', 1e-6, 'thr0', 1e-10, ...
    'NoiseLevel', 'none', 'eta', 1.01, 'RegParam0', 1, 'inSolver', 'gmres', ...
    'adaptConstr', 'sp', 'nonnegativity', 'off', 'verbosity', 'off',...
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
nn       = IRget(options, 'nonnegativity', [], 'fast');
inSolver = IRget(options, 'inSolver',      [], 'fast');

if strcmp(nn, 'on')
    nn = 1;
    options.adaptConstr = 'spnn';
else
    nn = 0;
end

if strcmp(inSolver, 'fgmres') && ~nn
    warning(['With options.inSolver = ''fgmres'' and options.nonnegativity = ''off'' ',...
        'it is more appropriate to use IRhybrid_fgmeres than to use IRirn.'])
end

% Call IRrestart with the specified options.
% Note that nonnegativity is specified via options.adaptConstr.
options = rmfield(options, 'nonnegativity');
[X, info] = IRrestart(A, b, K, options);