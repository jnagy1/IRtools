function val = PRget(options,name,value,flag)
%PRget Get options for IR Tools test problems
%
% val = PRget(options,'name',value,flag)
%
% val = PRget(options,'name') 
%   extracts the value of the named parameter from the options structure, 
%   returning an empty matrix if the parameter value is not specified in 
%   options.  It is sufficient to type only the leading characters that 
%   uniquely identify the parameter, and case is ignored for parameter names.  
%   [] is a valid options argument.
%
% val = PRget(options,'name',value) 
%   extracts the named parameter as above, but returns VALUE if the named 
%   parameter is not specified in options.
%
% val = PRget(options,'name',value,flag) 
%   with flag = 'fast' no error checking is performed.
%
% See also IRget, IRset, PRset

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

% Fast access with no error checking.
if (nargin == 4) && isequal(flag,'fast')
    val = PRgetfast(options,name,value);
    return
end

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3
    value = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('First argument must be an options structure created with PRset.');
end

if isempty(options)
    val = value;
    return;
end
allfields = {'trueImage';'PSF';'BlurLevel';'Frames';'BC';'CommitCrime';...
    'InterpMethod';'phantomImage';'CTtype';'sm';'angles';'p';'R';'d';...
    'span';'numCircles';'wavemodel';'s';'omega';'Tfinal';'Tsteps';...
    'material';'numData';'Tloglimits';'tauloglimits'};

Names = allfields;

name = deblank(name(:)'); % Force this to be a row vector.
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error(['Unrecognized property name ''%s''.  ' ...
        'See PRset for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        error('Ambiguous property name ''%s'' ', name);
    end
end

if any(strcmp(Names,Names{j,:}))
    val = options.(Names{j,:});
    if isempty(val)
        val = value;
    end
else
    val = value;
end

%------------------------------------------------------------------
function value = PRgetfast(options,name,defaultopt)
%PRGETFAST- Get PR OPTIONS parameter with no error checking.
%   VAL = PRGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.

if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct.
try
    value = options.(name);
catch
    value = [];
end

if isempty(value)
    value = defaultopt.(name);
end