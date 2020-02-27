function [X, info] = IRenrich(A, b, varargin)
% IRenrich Enriched CGLS method
%
% options  = IRenrich('defaults')
% [X,info] = IRenrich(A,b)
% [X,info] = IRenrich(A,b,K)
% [X,info] = IRenrich(A,b,options)
% [X,info] = IRenrich(A,b,K,options)
%
% IRenrich is a simplified driver for IRcgls with the enrichment option,
% giving the Enriched CGLS method (also known as Enriched CGNR).  This
% method finds the least squares solution in the enriched Krylov subspace
%   span{ A^Tb, A^TA A^Tb, (A^TA)^2 A^Tb, ..., (A^TA)^{k-1} A^Tb } + W ,
% where the subspace W is spanned by the columns of a user-defined
% enrichment matrix; the default is ones(N,1)/sqrt(N);
%
% It is also possible to solve the Tikhonov problem in standard from using
% the same enriched subspace.
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
%                   [ {'none'} | nonnegative scalar ]
%      eta        - safety factor for the discrepancy principle
%                   [ {1.01} | scalar greater than (and close to) 1 ]
%      NE_Rtol    - relative tolerance on the normal equation residual norm
%                   [ {1e-12} | positive integer ]
%      enrichment - matrix with a basis for the enrichment subspace
%                   [ {'ones'} | matrix ]
%      RegParam   - regularization parameter lambda, to be employed if
%                   IRenrich is used to solve the regularized problem
%                   (A'*A + lambda^2*L'*L)*x = A'*b;
%                   [ {0} | nonnegative scalar ]
%      RegMatrix  - regularization matrix L, used either as a
%                   regularization matrix to solve the regularized
%                   problem (A'*A + lambda^2*L'*L)*x = A'*b
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
%      its       - number of the last computed iteration
%      saved_iterations - iteration numbers of iterates stored in X 
%      StopFlag - a string that describes the stopping condition:
%                   * Reached maximum number of iterations
%                   * Residual tolerance satisfied (discrepancy principle) 
%                   * Normal equation residual tolerance satisfied
%      Rnrm      - relative residual norms at each iteration at each iteration
%      NE_Rnrm   - normal eqs relative residual norms at each iteration
%      Xnrm      - solution norms at each iteration
%      Enrm      - relative error norms (requires x_true) at each iteration
%      StopReg   - struct containing information about the solution that
%                  satisfies the stopping criterion.  Fields:
%                    It      : iteration where stopping criterion is satisfied
%                    X       : solution satisfying the stopping criterion
%                    Enrm    : the corresponding relative error (requires x_true)
%                    Rnrm    : enriched relative residual
%                    NE_Rnrm : enriched normal equations relative residual
%                    Xnrm    : solution norm
%      BestReg   - struct containing information about the solution that
%                  minimizes Enrm (requires x_true). Fields:
%                    It      : iteration where the minimum is attained
%                    X       : best solution
%                    Enrm    : best relative error
%                    Rnrm    : enriched relative residual for best solution
%                    NE_Rnrm : enriched normal equations relative residual
%                              for best solution
%                    Xnrm    : solution norm for best solution
%      StdCGLS   - struct containing information about the standard CGLS
%
% See also: IRcgls, IRget, IRset

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
    'enrichment','ones', 'stdCGLS_out','off', 'NoStop','off', ...
    'RegParam',0, 'RegMatrix','Identity', 'Reorth','on');
  
% If input is 'defaults,' return the default options in X.
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    X = defaultopt;
    return;
end

defaultopt.RegParam = 0;
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
    otherwise
        error('Too many input parameters')
end

if isempty(options)
    options = defaultopt;
end

% Call to IRcgls with the specified options.
options = IRset(defaultopt, options);
options = rmfield(options, 'MaxIter');
[X, info] = IRcgls(A, b, K, options);