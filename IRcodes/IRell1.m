function [X, info] = IRell1(A, b, varargin)
% IRell1 Least squares solver with 1-norm penalization term
%
% options  = IRell1('defaults')
% [X,info] = IRell1(A,b)
% [X,info] = IRell1(A,b,K)
% [X,info] = IRell1(A,b,options)
% [X,info] = IRell1(A,b,K,options)
%
% This penalized iterative method enforces sparsity on the solution via a
% 1-norm penalization term.  This function is a simplified driver for
% IRhybrid_fgmres.
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
%                   [ non-negative scalar | {'gcv'} | 'discrep' ]
%                   This also determines which stopping rule is used
%                   If 'gcv' is chosen, the iteration is stopped when
%                     the GCV function minimum stabilizes or increases 
%                     within a certain window of iterations (see 'stopGCV',
%                     'FlatTol' and 'MinTol').
%                   If 'discrep' is chosen, and NoiseLevel is provided,
%                     then the discrepancy principle is used as stopping
%                     criterion (see 'NoiseLevel' and 'eta').
%      stopGCV    - stopping criterion for the iterations when GCV is used
%                   [ GCVvalues | {'resflat'} ]
%      FlatTol    - tolerance for detecting flatness (stabilization)
%                   in the GCV function as a stopping criterion
%                   [ {10^-6} | non-negative scalar ]
%      MinTol     - window of iterations: if the GCV minimum continues
%                   to increase over this window, then the iterations are
%                   stopped:
%                   [ {3} | positive integer ]
%      NoiseLevel - norm of noise in rhs divided by norm of rhs
%                   (must be assigned if RegParam is 'discrep')
%                   [ {'none'} | nonnegative scalar ]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      RegParam0  - regularization parameter used in the first 
%                   projected problem (needed if RegParam is 'discrep')
%                   [ {1} | positive scalar ]
%      NoStop     - specifies whether the iterations should proceed
%                   after a stopping criterion has been satisfied
%                   [ 'on' | {'off'} ]
%      IterBar    - shows the progress of the outer iterations
%                   [ {'on'} | 'off' ]
% Note: the options structure can be created using the function IRset.
%
% Outputs:
%   X : computed solutions, stored column-wise (at the iterations listed in K)
%   info: structure with the following fields:
%      its          - number of the last computed iteration 
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag_in  - string that describes the inner stopping condition:
%                       * Stopping criterion of inner iterations is never
%                         satisfied
%                       * Stopping criterion is satisfied at least once
%                         during the inner iterations
%      StopFlag_out - string that describes the outer stopping condition:
%                       * Outer stopping criterion is never satisfied
%                       * Diagonal weighting matrix is numerically zero
%                       * Solution stabilizes
%                       * Transformed solution stabilizes
%                       * Regularization parameter stabilizes
%      Rnrm    - relative residual norms at each iteration
%      Xnrm    - solution norms at each iteration
%      Enrm    - relative error norms (requires x_true) at each iteration
%      StopReg - struct containing information about the solution that
%                satisfies the stopping criterion.  Fields:
%                  It   : iteration where the stopping criterion is satisfied
%                  X    : solution satisfying the stopping criterion
%                  Enrm : the corresponding relative error (requires x_true)
%      BestReg - struct containing information about the solution that
%                minimizes Enrm (requires x_true).  Fields:
%                  It   : iteration where the minimum is attained
%                  X    : best solution
%                  Enrm : best relative error
%
% See also: IRhtv, IRhybrid_fgmres, IRirn, IRget, IRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Set default values for options.
defaultopt = struct('x0', 'none', 'MaxIter', 100 , 'RegParam', 'gcv', ...
    'stopGCV', 'resflat', 'resflatTol', 0.05, 'GCVflatTol', 10^-6, ...
    'GCVminTol', 3, 'x_true', 'none', 'IterBar', 'on', 'NoStop', 'off', ...
    'thr0', 1e-10, 'NoiseLevel', 'none', 'eta', 1.01, 'RegParam0', 1);
  
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

% Call IRhybrid_fgmres with the specific options.
options = IRset(defaultopt, options);
if isempty(K)
    [X, info] = IRhybrid_fgmres(A, b, options);
else
    [X, info] = IRhybrid_fgmres(A, b, K, options);
end