function vout = Ptransp_solve(L, vin)
%Ptransp_solve  Auxi,ilary function for IR Tools
%
% vout = Ptransp_solve
%
% This computes vout = L'\vin.  If L is a function handle, then we use the
% user-supplied function, which is passed via options as 'RegMatrix'.
% Otherwise, we use the backslash operator.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

if isa(L, 'function_handle')
    transp_flag = 'transp';
    vout = L(vin, transp_flag);
else
    vout = (L')\vin;
end