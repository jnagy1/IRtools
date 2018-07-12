function options = PRset(varargin)
%PRset Set options for IR Tools test problems
%
% options = PRset('param1',value1,'param2',value2,...)
%
% Create/alter options structure for test problems in IT Tools, in which
% the named parameters have the specified values.  Any unspecified
% parameters are set to [] (they indicate to use the default value for that 
% parameter when passed to an test problem). It is sufficient to type only
% the leading characters that uniquely identify the parameter, and case is 
% ignored for parameter names.
% NOTE: For values that are strings, the complete string is required.
%
% options = PRset(OLDOPTS,'param1',value1,...) creates a copy of OLDOPTS
% with the named parameters altered with the specified values.
%
% options = PRset(OLDOPTS,NEWOPTS) combines an existing options structure
% OLDOPTS with a new structure NEWOPTS.  Any parameters in NEWOPTS with
% non-empty values overwrite the corresponding old parameters in OLDOPTS.
%
% PRset with no input arguments and no output arguments displays all
% parameter names and their possible values, with defaults shown in {}.
%
% options = PRset(with no input arguments) creates an options structure
% where all the fields are set to [].
%
% options = PRset('PRproblem') creates an options structure with all the
% parameter names and default values relevant to an IR method. That is,
%           PRset('PRproblem')
% or
%           PRset(@PRproblem)
% returns an options structure containing all the parameter names and
% default values relevant to a test problem.
%
% See also IRget, IRset, PRget

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Print out possible values of properties. 
if (nargin == 0) && (nargout == 0)
  fprintf('PRblur       trueImage    [ array | ''pattern1'' | ''pattern2'' | ''sppattern'' | ''ppower'' |\n')
  fprintf('                            ''ppower'' | ''dot2'' | ''dotk'' | ''satellite'' | {''hst''} ]\n')
  fprintf('PRblur       PSF          [ array | {''gauss''} | ''defocus'' | ''speckle'' | ''shake'' |\n')
  fprintf('                            ''motion'' | ''rotation'' ]\n')
  fprintf('PRblur       BlurLevel    [ ''mild'' | {''medium''} | ''severe'' ]\n')
  fprintf('PRblur       BC           [ ''zero'' | ''periodic'' | {''reflective''} ]\n')
  fprintf('PRblur       CommitCrime  [ {''off''} | ''on'' ]\n')
  fprintf('PRdiffusion  Tfinal       [ positive scalar | {0.01} ]\n')
  fprintf('PRdiffusion  Tsteps       [ positive integer | {100} ]\n')
  fprintf('PRinvinterp2 InterpMethod [ ''nearest'' | {''linear''} | ''cubic'' | ''spline'' ]\n')
  fprintf('PRnmr        numData      [ {''double'' | m | [m1, m2] ]\n')
  fprintf('PRnmr        material     [ {''carbonate''} | ''methane'' | ''organic'' | ''hydroxyl'' ]\n')
  fprintf('PRnmr        numData      [ {''double'' | m | [m1, m2] ]\n')
  fprintf('PRnmr        Tloglimits   [ {[-4, 1]} | vector with two real components ]\n')
  fprintf('PRnmr        tauloglimits [ {[-4, 1]} | vector with two real components ]\n')
  fprintf('PRseismic    phantomImage [ {''tectonic''} | ''smooth'' | ''binary'' | ''threephases'' |\n')
  fprintf('                            ''threephasessmooth'' | '' fourphases'' | ''grains'' | ''ppower'' ]\n')
  fprintf('PRseismic    wavemodel    [ {''ray''} | ''fresnel'' ]\n')
  fprintf('PRseismic    s            [ positive integer | {n} ]\n')
  fprintf('PRseismic    p            [ positive integer | {2*n} ]\n')
  fprintf('PRseismic    omega        [ positive scalar | {10} ]\n')
  fprintf('PRseismic    sm           [ {''true''} | false ]\n')
  fprintf('PRspherical  phantomImage [ {''shepplogan''} | ''smooth'' | ''binary'' | ''threephases'' |\n')
  fprintf('                            ''threephasessmooth'' | '' fourphases'' | ''grains'' | ''ppower'' ]\n')
  fprintf('PRspherical  angles       [ vector with positive scalars | {linspace(360/n,360,n)} ]\n');
  fprintf('PRspherical  numCircles   [ positive integer | {round(sqrt(2)*n)} ]\n');
  fprintf('PRseismic    sm           [ {''true'' | false ]\n');
  fprintf('PRtomo       phantomImage [ {''shepplogan''} | ''smooth'' | ''binary'' | ''threephases'' | |\n')
  fprintf('                            ''threephasessmooth'' | '' fourphases'' | ''grains'' | ''ppower'' ]\n');
  fprintf('PRtomo       CTtype       [ {''parallel''} | ''fancurved'' ]\n');
  fprintf('PRtomo       sm           [ {logical true} | logical false]\n');
  fprintf('PRtomo       angles       [ vector of positive scalars | |\n')
  fprintf('                            {0:1:179 (for ''parallel''),0:2:258 (for ''fancurved'')} ]\n');
  fprintf('PRtomo       p            [ positive integer | {round(sqrt(2)*n)} ]\n');
  fprintf('PRtomo       d            [ positive scalar | {d-1} ]\n');
  fprintf('PRtomo       R            [ positive scalar | {2} ]\n');
  fprintf('PRtomo       span         [ positive scalar | {see documentation} ]\n');
  return
end

% Create a struct of all the fields.
allfields = {'trueImage';'PSF';'BlurLevel';'Frames';'BC';'CommitCrime';...
    'InterpMethod';'phantomImage';'CTtype';'sm';'angles';'p';'R';'d';...
    'span';'numCircles';'wavemodel';'s';'omega';'Tfinal';'Tsteps';...
    'material';'numData';'Tloglimits';'tauloglimits'};
  
% Create cell array.
structinput = cell(2,length(allfields));
% Fields go in first row.
structinput(1,:) = allfields';
% []'s go in second row.
structinput(2,:) = {[]};
% Turn it into correctly ordered comma separated list and call struct.
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
     error('Undefined function.  Please use PRset(''funcname'') or PRset(@funcname), where funcname is a PR problem')
  end
  try 
    optionsfcn = feval(varargin{1},'defaults');
  catch
    error('PRset ONLY works with an PR method.  Please use PRset(@PR).')
  end
  % The defaults from IR functions don't include all the field typically,
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
  if ischar(arg)                        % arg is an option name.
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument.
    if ~isa(arg,'struct')
      error(['Expected argument %d to be a string parameter name ' ...
          'or an options structure \n created with PRset.'], i);
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

expectval = 0;                % Start expecting a name, not a value.
while i <= numberargs
  arg = varargin{i};

  if ~expectval
    if ~ischar(arg)
      error('Expected argument %d to be a string parameter name.', i);
    end

    lowArg = lower(arg);
    %j = strmatch(lowArg,names);
    j = find(strcmp(lowArg, names));
    if isempty(j)                       % If no matches.
      error('Invalid parameter name ''%s'' ', arg);
    elseif length(j) > 1                % If more than one match.
      % Check for any exact matches (in case any names are subsets of others)
      %k = strmatch(lowArg,names);
      k = find(strcmp(lowArg, names));
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


%----SUBFUNCTION---------------------------------------------
function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
  return
end

% See if 'field' is a valid field.
%validfield = true;
switch field
  case {'trueImage'} % numeric array, various char strings
    [validvalue, errmsg] = trueImageType(field,value);
  case {'PSF'} % numeric array, various char strings
    [validvalue, errmsg] = psftype(field,value);
  case {'BlurLevel'} %  mild, medium of severe
    [validvalue, errmsg] = BlurLevelType(field,value);
  case {'Frames'} % positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case {'BC'} %  zero, reflective, neumann, periodic
    [validvalue, errmsg] = BCtype(field,value);
  case {'CommitCrime'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case{'InterpMethod'} % nearest, linear, cubic, spline 
    [validvalue, errmsg] = InterpMethodType(field,value);
  case {'phantomImage'} % numeric array, various char strings
    [validvalue, errmsg] = phantomType(field,value);
   case {'CTtype'} %  parallel or fancurved
    [validvalue, errmsg] = CTtypeType(field,value);  %%%%%%%
  case {'wavemodel'} %  ray or fresnel %%%%%%%
    [validvalue, errmsg] = wavemodelType(field,value);
  case {'sm'} % logical true or false 
    [validvalue, errmsg] = TrueFalseType(field,value);
  case {'p'} % real positive scalar
    [validvalue, errmsg] = positiveScalarNaN(field,value);
  case {'R'} % real positive scalar
    [validvalue, errmsg] = positiveScalar(field,value);
  case {'d'} % real positive scalar
    [validvalue, errmsg] = positiveScalarNaN(field,value);
  case {'span'} % real positive scalar
    [validvalue, errmsg] = positiveScalarNaN(field,value);
  case {'angles'} % numeric vector
    [validvalue, errmsg] = RealVectorType(field,value);
  case {'numCircles'} % positive integer
    [validvalue, errmsg] = positiveScalarNaN(field,value);
  case {'s'} % real positive scalar
    [validvalue, errmsg] = positiveScalar(field,value);
  case {'omega'} % positive integer
    [validvalue, errmsg] = positiveScalar(field,value);
  case {'Tfinal'} % positive integer
    [validvalue, errmsg] = positiveScalar(field,value);
  case {'Tsteps'} % positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case {'material'} % numeric array, various char strings
    [validvalue, errmsg] = materialType(field,value);
  case {'numData'} % numeric array, various char strings
    [validvalue, errmsg] = numDataType(field,value);
  case {'Tloglimits'} % numeric vector
    [validvalue, errmsg] = VectLimitType(field,value);
  case {'tauloglimits'} % numeric vector
    [validvalue, errmsg] = VectLimitType(field,value);
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

%-----------------------------------------------------------------------

function [valid, errmsg] = trueImageType(field,value)
% Either a numeric array or char string
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'pattern1';'pattern2';'sppattern';'ppower';'smooth';'dot2';'dotk';'satellite';'hst'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''pattern1'' or ''pattern2'' or ''ppower'' or ''dot2'' or ''dotk'' or ''satellite'' or ''hst''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = psftype(field,value)
% Either a numeric array or a char string
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'gauss';'defocus';'rotation';'speckle';'motion';'shake'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''gauss'' or ''defocus'' or ''speckle'' or ''rotation'' or ''motion'' or ''shake''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = BlurLevelType(field,value)
% One of these strings: mild, medium, severe
valid =  ischar(value) && any(strcmp(value,{'mild';'medium';'severe'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''mild'' or ''medium'' or ''severe''.',field);
else
  errmsg = '';
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

function [valid, errmsg] = BCtype(field,value)
% One of the given strings. Note that reflective = neumann = reflexive
valid =  ischar(value) && any(strcmp(value,{'zero';'reflective';'neumann';'reflexive';'periodic';'crime'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''zero'' or ''reflective'' (or neuman or reflexive) or ''periodic''.',field);
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

%------------------------------------------------------------------------

function [valid, errmsg] = InterpMethodType(field,value)
% One of the given strings: nearest, linear, cubic, spline
valid =  ischar(value) && any(strcmp(value,{'nearest';'linear';'cubic';'spline'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''nearest'' or ''linear'' or ''cubic'' or ''spline''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = phantomType(field,value)
% Either a numeric array or char string
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'shepplogan';'tectonic';'smooth';'binary';'threephases';'threephasessmooth';'fourphases';'grains';'ppower'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''shepplogan'', ''tectonic'', ''smooth'', ''binary'', ''threephases'', ''threephasessmooth'', ''fourphases'', ''grains'' or ''ppower''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = RealVectorType(field,value)
% Any vector with real numeric values
valid =  isnumeric(value) && isreal(value) && numel(value)==max(size(value));

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be real values.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = positiveScalar(field,value)
% Any real non-negative scalar
valid =  isreal(value) && isscalar(value) && (value > 0);

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = CTtypeType(field,value)
% One of the given strings. 
valid =  ischar(value) && any(strcmp(value,{'parallel';'fancurved'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''parallel'' or ''fancurved''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = TrueFalseType(field,value)
% One of these strings: on, off
valid =  isa(value,'logical');
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be logical true or false.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = positiveScalarNaN(field,value)
% Any real non-negative scalar or NaN
valid =  (isreal(value) && isscalar(value) && (value > 0)) || isnan(value);

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real positive scalar.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = wavemodelType(field,value)
% One of the given strings. 
valid =  ischar(value) && any(strcmp(value,{'ray';'fresnel'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''ray'' or ''fresnel''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = materialType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value, {'carbonate', 'methane', 'organic', 'water'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''carbonate'', ''methane'', ''organic'', ''water''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = numDataType(field,value)
% One of these strings: mild, medium, severe
valid =  (ischar(value) && strcmp(value,'double')) ||...
         (isreal(value) && (all(value > 0)) && all(value == floor(value))...
          && (numel(value) == 1 || numel(value) == 2));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''mild'' or ''medium'' or ''severe''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = VectLimitType(field,value)
% Any vector with real numeric values
valid =  isnumeric(value) && isreal(value) && numel(value)==2;

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be two real values.',field);
else
  errmsg = '';
end