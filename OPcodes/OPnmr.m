function w = OPnmr(u,A1,A2,tflag)
% OPnmr The forward computation and its adjoint for PRnmr
%
% w = OPnmr(u,A1,A2,tflag)
%
% Performs either the forward computation, its ajoint, or returns the
% problem dimensions.
%
% Input: u     - the vector to be operated upon
%        A1,A2 - matrices for 1D kernels with respect to the (tau1, T1)
%                and (tau2, T2) variables, respectively
%        tflag - string that determines the computation:
%                'notransp', 'transp' or 'size'

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD Licence. A separate license file should be provided as part 
% of the package.

[m1,n1] = size(A1);
[m2,n2] = size(A2);
if strcmpi(tflag,'size')
    w(1) = m1*m2;
    w(2) = n1*n2;
elseif strcmpi(tflag,'notransp')
    u = reshape(u, n1, n2);
    w = A1*u*A2';
    w = w(:);
else
    u = reshape(u, m1, m2);
    w = A1'*u*A2;
    w = w(:);
end