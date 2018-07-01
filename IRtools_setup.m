function IRtools_setup
%IRtools_setup Set up search paths to IR Tools
%
%  Run this function to setup IR Tools. 
%
%  In most cases, this script should only need to be run once, and then the 
%  path will be permanently saved. Note the following:
%    - It forces all user paths in the current session to set permanently.
%    - On systems with shared MATLAB licenses, a permanent save may need to
%      be done manually. In this case, see SAVEPATH for more information.
%  
%  If you want to remove the paths associated with IR Tools, use RMPATH.
%
% See also: addpath, rmpath, savepath

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

% Add IR Tools to the search path.

addpath(genpath(fileparts(mfilename('fullpath'))));

% Check if AIR Tools II is installed, and issue a warning if it is not or
% if an old version is installed.
%
% Note that recent versions of MATLAB recommend to use ~contains in the next
% if statement, but this does not work for some older versions of MATLAB.
% So we use isempty.
if isempty(strfind(lower(path),'airtools'))
    warning('Could not find AIR Tools II; some functionality may be limited.')
    type('Extra/INFO_AIRToolsII_not_installed.txt')
else
    if exist('ARTdemo.m','file')
        if exist('AIRToolsII_setup.m','file')
            AIRtools_old_path = which('ARTdemo.m');
            AIRtools_old_path = AIRtools_old_path(1:end-9);
            rmpath(AIRtools_old_path);
            warning('Removing old AIR Tools from path, and using AIR Tools II')
        else
            warning('An old version of AIR Tools is installed; this version will not work with IR Tools.')
            type('INFO_old_AIRTools_installed.txt')
        end
    end
end
status = savepath;
if status == 1
    warning('IR Tools was added to the MATLAB search path for the current session only. Adding it permanently failed, probably due to a write permission issue. It is possible to manually add IR Tools permanently to your search path, but may require consulting your system administrator. Alternatively, you can re-run this installation function in each new MATLAB session where you want to use IR Tools.')
end
