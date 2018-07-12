function options = IRset(varargin)
%IRset Set options for IR Tools functions
%
% options = IRset('param1',value1,'param2',value2,...)
%
% Create/alter options structure for methods in IR Tools, in which the
% named parameters have the specified values.  Any unspecified parameters
% are set to [] (they indicate to use the default value for that parameter
% when passed to an IR function). It is sufficient to type only the
% leading characters that uniquely identify the parameter, and  case is
% ignored for parameter names.
% NOTE: For values that are strings, the complete string is required.
%
% options = IRset(OLDOPTS,'param1',value1,...) creates a copy of OLDOPTS
%   with the named parameters altered with the specified values.
%
% options = IRset(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new structure NEWOPTS.  Any parameters in NEWOPTS with
%   non-empty values overwrite the corresponding old parameters in OLDOPTS.
%
% options = IRset('IRmethod') creates an options structure with all the
% parameter names and default values relevant to an IR method. That is,
%           IRset('IRmethod')
% or
%           IRset(@IRmethod)
% returns an options structure containing all the parameter names and
% default values relevant to an IR method.
%
% Examples:
%   To create options with the default options for an IR nethod
%       options = IRset('IRmethod');
%   To change the maximum iterations to 150 in options
%       options = IRset(options,'MaxIter',150);
%
% See also: IRget, PRget, PRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Create a struct of all the fields. Note that if new options are added
% then IRoptions_fields.m should be edited to include the field name in
% the list created in IRoptions_fields.m.
%
allfields = IRoption_fields;
  
% Create cell array.
structinput = cell(2,length(allfields));
% Fields go in first row.
structinput(1,:) = allfields';
% []'s go in second row.
structinput(2,:) = {[]};
% Turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % We might change this value, so assign it.

% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
  if ischar(varargin{1})
    funcname = lower(varargin{1});
  elseif isa(varargin{1},'function_handle')
    funcname = func2str(varargin{1});
  end
  if ~exist(funcname, 'file')
     error('Undefined function.  Please use IRset(''funcname'') or IRset(@funcname), where funcname is a IR method')
  end
  try 
    optionsfcn = feval(varargin{1},'defaults');
  catch
    error('IRset ONLY works with an IR method.  Please use IRset(@IR).')
  end
  % The defaults from IR functions don't include all the fields typically,
  % so run the rest of IRset as if called with IRset(options,optionsfcn)
  % to get all the fields.
  varargin{1} = options;
  varargin{2} = optionsfcn;
  numberargs = 2;
end

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name.
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument.
    if ~isa(arg,'struct')
      error(['Expected argument %d to be a string parameter name ' ...
             'or an options structure \n created with IRset.'], i);
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),Names{j,:}))
        val = arg.(Names{j,:});
      else
        val = [];
      end
      if ~isempty(val)
        if ischar(val)
          val = lower(deblank(val));
        end
        checkfield(Names{j,:},val)
        options.(Names{j,:}) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end

expectval = 0;                    % Start expecting a name, not a value.
while i <= numberargs
  arg = varargin{i};

  if ~expectval
    if ~ischar(arg)
      error('Expected argument %d to be a string parameter name.', i);
    end

    lowArg = lower(arg);
    j = find(strcmp(lowArg,names));
    if isempty(j)                       % If no matches.
      error('Invalid parameter name ''%s'' ', arg);
    elseif length(j) > 1                % If more than one match.
      % Check for any exact matches (in case any names are subsets of others).
      k = find(strcmp(lowArg,names));
      if length(k) == 1
        j = k;
      else
        error('Ambiguous parameter name ''%s'' ', arg);
      end
    end
    expectval = 1;                      % We expect a value next.

  else
    if ischar(arg)
      arg = lower(deblank(arg));
    end
    checkfield(Names{j,:},arg);
    options.(Names{j,:}) = arg;
    expectval = 0;
  end
  i = i + 1;
end

if expectval
  error( 'Expected value for parameter ''%s''.', arg);
end


%----SUBFUNCTIONS---------------------------------------------

function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 

% Empty matrix is always valid.
if isempty(value)
  return
end

% See if 'field' is a valid field.
%validfield = true;
switch field
  case {'MaxIter'} % real positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case {'MaxIterIn'} % real positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case{'TotIterMax'}
    [validvalue, errmsg] = PosInteger2(field, value);
  case{'ktotcount'}
    [validvalue, errmsg] = nonNegInteger2(field, value);
  case {'MaxIterOut'} % real positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case {'x0'} % numeric array, none
    [validvalue, errmsg] = x0type(field,value);
  case {'IterBar'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'x_true'} % numeric array, off
    [validvalue, errmsg] = x_truetype(field,value);
  case {'NoStop'} %  off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'NoStopIn'} %  off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'NoStopOut'} %  off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'restart'} %  off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'MatUpdate'}
    [validvalue, errmsg] = onOffType(field,value);
  case {'stopGCV'}
    [validvalue, errmsg] = stopGCVtype(field,value);
  case {'stopCrit'}
    [validvalue, errmsg] = stopCritype(field,value);
  case {'stopOut'}
    [validvalue, errmsg] = stopOutype(field,value);
  case {'tollalpha'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'NE_Rtol'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'thr0'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'stabOut'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'RegMatrix'} % character string
    [validvalue, errmsg] = tikFunType(field,value);
  case {'enrichment'} % character string
    [validvalue, errmsg] = enrichType(field,value);
  case {'stdCGLS_out'}
    [validvalue, errmsg] = onOffType(field,value);
  case {'transf'} % character string
    [validvalue, errmsg] = spTransf(field,value);
  case {'Reorth'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'stoprule'} % character string
    [validvalue, errmsg] = stopruleType(field,value);
  case {'xMin'} % real scalar
    [validvalue, errmsg] = realScalar(field,value);
  case {'xMax'} % real scalar
    [validvalue, errmsg] = realScalar(field,value);
  case {'xEnergy'} % real positive scalar or 'none'
    [validvalue, errmsg] = positiveScalar2(field,value);
  case {'taudelta'} % real non-negative scalar or 'none'
    [validvalue, errmsg] = nonNegScalar2(field,value);
  case {'nonneg'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'verbosity'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'Ubound'} % real positive scalar or 'none'
    [validvalue, errmsg] = positiveScalar2(field,value);
  case {'relaxParam'}% real non-negative scalar
    [validvalue, errmsg] = nonNegScalar2(field,value);
  case {'RegParam'} % 
    [validvalue, errmsg] = RegPartype(field,value);
  case {'GCVweight'} % real non-negative scalar, adapt
    [validvalue, errmsg] = WeightType(field,value);
  case {'GCVflatTol'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'tolX'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'GCVminTol'}% real positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case {'resflatTol'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'LSQRtols'}% real non-negative vector of length 2
    [validvalue, errmsg] = nonNeg2vector(field,value);
  case {'RegParam0'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'NoiseLevel'}% real non-negative scalar
    [validvalue, errmsg] = NoiseLevelType(field,value);
  case {'eta'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'BestEnrm'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'trunc'}
    [validvalue, errmsg] = truncation(field,value);
  case {'noiseType'}
    [validvalue, errmsg] = noiseType(field,value);
  case {'noiseParam'}
    [validvalue, errmsg] = PosVector2(field,value);
  case {'inSolver'} % character string
    [validvalue, errmsg] = inSolverType(field,value);
  case {'adaptConstr'} % character string
    [validvalue, errmsg] = adaptConstrType(field,value);
  case {'Ktot'} % vector with all real positive integers
    [validvalue, errmsg] = PosIntegerVector(field,value);
  case {'nonnegativity'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'DecompOut'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case{'sirt_method'}
    [validvalue, errmsg] = sirtType(field,value);  
  otherwise
    %validfield = false;  
    validvalue = false;
    errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
end

if validvalue 
    return;
else
  error(errmsg)
end

%------------------------------------------------------------------------

function [valid, errmsg] = PosIntegerVector(field,value)
% Any vector with all postive real integers
valid =  (isreal(value) && (all(value > 0)) && all(value == floor(value))) | (ischar(value) && any(strcmp(value,{'none'}))) ;
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be positive real integers.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = PosVector2(field,value)
% Any vector with all postive real integers
valid =  (isreal(value) && (all(value > 0))) && (length(value) == 2) || (ischar(value) && any(strcmp(value,{'none'}))) ;
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be positive real vector of length 2.',field);
else
    if value(1)>value(2)
        warning('For the methods to be meaningful, the Poisson parameter (%s(2)) should be bigger than the Gaussian std deviation (%s(1)).', field, field);
        errmsg = '';
    else
        errmsg = '';
    end
end

%------------------------------------------------------------------------

function [valid, errmsg] = PosInteger(field,value)
% Any positive real integer
valid =  isreal(value) && isscalar(value) && (value > 0) && value == floor(value) ;
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive real integer.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = PosInteger2(field,value)
% Any positive real integer
valid =  (isreal(value) && isscalar(value) && (value > 0) && value == floor(value)) || ... 
    (ischar(value) && strcmp(value,{'none'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive real integer, or ''none''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = nonNegInteger2(field,value)
% Any positive real integer
valid =  (isreal(value) && isscalar(value) && (value >= 0) && value == floor(value)) || ... 
    (ischar(value) && strcmp(value,{'none'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive real integer, or ''none''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = x_truetype(field,value)
% Either a numeric array or off
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'none'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''none''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = x0type(field,value)
% Either a numeric array or none
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'none'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''none''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = nonNegscalar(field,value)
% Any real non-negative scalar
valid =  isreal(value) && isscalar(value) && (value >= 0);
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = realScalar(field,value)
% Any real  scalar
valid =  (isreal(value) && isscalar(value)) || (ischar(value) && strcmpi(value, 'none'));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = tikFunType(field,value)
% One of these strings: none, Identity, Laplacian
valid =  (ischar(value) && any(strcmpi(value,{'none';'Identity';'Laplacian1D';'Laplacian2D';'tv'}))) || ...
         (~ischar(value) && ismatrix(value));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'' or ''Identity'' or ''Laplacian1D'' or ''Laplacian2D'' or ''tv'' or a matrix.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = enrichType(field,value)
% One of these strings: none, Identity, Laplacian
valid =  (ischar(value) && (strcmpi(value,'none')||strcmpi(value,'ones')))|| ...
         (~ischar(value) && ismatrix(value));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'' or ''ones'' or a matrix of coherent size.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = spTransf(field,value)
% One of these strings: none, Identity, Laplacian
valid =  ischar(value) && any(strcmpi(value,{'none';'dct'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'' or ''dct''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = stopruleType(field,value)
% One of these strings: none, Identity, Laplacian
valid =  ischar(value) && any(strcmpi(value,{'none';'NCP';'DP'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'' or ''NCP'' or ''DP''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = nonNegScalar2(field,value)
% Either a numeric array or none
valid =  (isreal(value) && isscalar(value) && value >= 0) || (ischar(value) && any(strcmp(value,{'none'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be non-negative scalar or ''none''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = positiveScalar2(field,value)
% Either a numeric array or none
valid =  (isreal(value) && isscalar(value) && value > 0) || (ischar(value) && any(strcmp(value,{'none'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be positive scalar or ''none''.',field);
else
  errmsg = '';
end%--------------------------------------------------------------------------

function [valid, errmsg] = RegPartype(field,value)
% One of these: real nonnegative scalar, GCV, WGCV, optimal
valid =  (isreal(value) && isscalar(value) && (value >= 0)) | (ischar(value) && any(strcmpi(value,{'gcv','wgcv','modgcv','optimal','discrep','off'})));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''GCV'' or ''WGCV'' or ''discrep'' or ''optimal''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = stopGCVtype(field,value)
% One of these: resflat, GCVflat
valid =  (ischar(value) && any(strcmpi(value,{'resflat', 'GCVvalues'})));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''resflat'' or ''GCVvalues''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = stopCritype(field,value)
% One of these: resflat, GCVflat
valid =  (ischar(value) && any(strcmpi(value,{'none', 'resflat', 'discrep'}))); % 'outer'
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''resflat'' or ''discrep''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = stopOutype(field,value)
% One of these: xstab, Lxstab, regPstab
valid =  (ischar(value) && any(strcmpi(value,{'xstab', 'Lxstab', 'regPstab'})));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''xstab'' or ''Lxstab'' or ''regPstab''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = WeightType(field,value)
% One of these: real non-negative scalar, adapt
valid =  (isreal(value) && isscalar(value) && (value >= 0) )| (ischar(value) && any(strcmp(value,{'adapt'})));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''adapt''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

%function [valid, errmsg] = InSolvetype(field,value)
% % One of these strings: tsvd, tikhonov, none
%valid =  ischar(value) && any(strcmp(value,{'tsvd','tikhonov','none'}));
%if ~valid
% errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''tsvd'', ''tikhonov'' or ''none''.',field);
%else
%  errmsg = '';
%end

%------------------------------------------------------------------------

function [valid, errmsg] = nonNeg2vector(field,value)
% Any real non-negative scalar
valid =  isreal(value) && isvector(value) && sum(value>=0)==length(value) && (length(value)==2);
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative vector of length 2.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = NoiseLevelType(field,value)
% One of these: real nonnegative scalar, GCV, WGCV, optimal
valid =  (isreal(value) && isscalar(value) && (value >= 0)) |...
         (ischar(value));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: a nonnegative scalar noise level should be provided.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = truncation(field,value)
% One of these: real nonnegative scalar, GCV, WGCV, optimal
valid =  (isreal(value) && isscalar(value) && (value > 0) && value == floor(value)) | (ischar(value) && any(strcmpi(value,{'full', 'none'})));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive integer or ''none'' or ''full''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = noiseType(field,value)
% One of these strings: none, Identity, Laplacian
valid =  ischar(value) && any(strcmpi(value,{'gauss';'poisson'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''gauss'' or ''poisson''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = inSolverType(field,value)
% One of these strings: gmres, lsqr, fgmres
valid =  (ischar(value) && any(strcmpi(value,{'gmres';'fgmres';'lsqr';'cgls';'rrgmres'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS inSolver %s: must be ''gmres'' or ''fgmres'' or ''lasq'' or ''cgls''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = adaptConstrType(field,value)
% One of these strings: tv, nn
valid =  (ischar(value) && any(strcmpi(value,{'tv';'nn';'tvnn';'box';'energy';'project';'spnn';'sp'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS inSolver %s: must be ''tv'' or ''nn'' or ''tvnn'' or ''box'' or ''spnn''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = sirtType(field,value)
% One of these strings: cav, cimmino, drop, landweber, sart
valid =  (ischar(value) && any(strcmpi(value,{'cav';'cimmino';'drop';'landweber';'sart'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''cav'' or ''cimmino'' or ''drop'' or ''landweber'' or ''sart''.',field);
else
  errmsg = '';
end