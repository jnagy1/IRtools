function vout = Atransp_times_vec(A, vin)
%Atransp_times_vect  Auxiliary function for IR Tools
%
% This computes vout = A'*vin.  If A is a function handle, then we use the
% user-supplied function, which is passed as A. Otherwise, we use the
% mtimes * operator.
%
% See also: A_times_vec

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

if isa(A, 'function_handle')
    transp_flag = 'transp';
    vout = A(vin, transp_flag);
else
    vout = A'*vin;
end